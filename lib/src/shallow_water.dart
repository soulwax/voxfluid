import 'dart:math' as math;
import 'dart:typed_data';

/// Eulerian 2D shallow-water virtual pipe simulator conserving fluid volume across uneven terrain.
final class ShallowWaterGrid {
  final int width;
  final int height;
  final double cellSize;
  double gravity;
  double pipeCrossSection;
  double damping;
  double evaporationRate;

  final Float64List _terrain;
  final Float64List _depth;
  final Float64List _nextDepth;
  final Float64List _fluxL;
  final Float64List _fluxR;
  final Float64List _fluxU;
  final Float64List _fluxD;
  final Float64List _vx;
  final Float64List _vz;

  ShallowWaterGrid(
    this.width,
    this.height, {
    this.cellSize = 1.0,
    this.gravity = 9.81,
    this.pipeCrossSection = 1.0,
    this.damping = 0.992,
    this.evaporationRate = 0.0,
  })  : _terrain = Float64List(width * height),
        _depth = Float64List(width * height),
        _nextDepth = Float64List(width * height),
        _fluxL = Float64List(width * height),
        _fluxR = Float64List(width * height),
        _fluxU = Float64List(width * height),
        _fluxD = Float64List(width * height),
        _vx = Float64List(width * height),
        _vz = Float64List(width * height);

  int _idx(int x, int z) => z * width + x;
  bool inBounds(int x, int z) => x >= 0 && x < width && z >= 0 && z < height;

  double terrainAt(int x, int z) => inBounds(x, z) ? _terrain[_idx(x, z)] : 0.0;

  void setTerrainHeight(int x, int z, double h) {
    if (inBounds(x, z)) _terrain[_idx(x, z)] = h;
  }

  void fillTerrainFromSampler(
    double Function(double worldX, double worldZ) sampler, {
    double originX = 0.0,
    double originZ = 0.0,
  }) {
    for (var z = 0; z < height; z++) {
      for (var x = 0; x < width; x++) {
        _terrain[_idx(x, z)] = sampler(originX + x * cellSize, originZ + z * cellSize);
      }
    }
  }

  double waterDepthAt(int x, int z) => inBounds(x, z) ? _depth[_idx(x, z)] : 0.0;

  double surfaceHeightAt(int x, int z) =>
      inBounds(x, z) ? (_terrain[_idx(x, z)] + _depth[_idx(x, z)]) : 0.0;

  ({double vx, double vz}) velocityAt(int x, int z) {
    if (!inBounds(x, z)) return (vx: 0.0, vz: 0.0);
    final i = _idx(x, z);
    return (vx: _vx[i], vz: _vz[i]);
  }

  void addWaterDepth(int x, int z, double depthMeters) {
    if (!inBounds(x, z) || depthMeters <= 0.0) return;
    _depth[_idx(x, z)] += depthMeters;
  }

  void addWaterVolume(int x, int z, double cubicMeters) {
    if (!inBounds(x, z) || cubicMeters <= 0.0) return;
    final area = math.max(1e-6, cellSize * cellSize);
    _depth[_idx(x, z)] += cubicMeters / area;
  }

  void clearWater() {
    _depth.fillRange(0, _depth.length, 0.0);
    _nextDepth.fillRange(0, _nextDepth.length, 0.0);
    _fluxL.fillRange(0, _fluxL.length, 0.0);
    _fluxR.fillRange(0, _fluxR.length, 0.0);
    _fluxU.fillRange(0, _fluxU.length, 0.0);
    _fluxD.fillRange(0, _fluxD.length, 0.0);
    _vx.fillRange(0, _vx.length, 0.0);
    _vz.fillRange(0, _vz.length, 0.0);
  }

  double get totalWaterVolume {
    final area = cellSize * cellSize;
    var sum = 0.0;
    for (var i = 0; i < _depth.length; i++) {
      sum += _depth[i] * area;
    }
    return sum;
  }

  int get floodedCellCount {
    var c = 0;
    for (var i = 0; i < _depth.length; i++) {
      if (_depth[i] > 0.005) c++;
    }
    return c;
  }

  double get maxWaterDepth {
    var m = 0.0;
    for (var i = 0; i < _depth.length; i++) {
      if (_depth[i] > m) m = _depth[i];
    }
    return m;
  }

  double sampleSurfaceHeightNormalized(double nx, double nz) {
    final gx = nx.clamp(0.0, 1.0) * math.max(0, width - 1);
    final gz = nz.clamp(0.0, 1.0) * math.max(0, height - 1);
    final x0 = gx.floor().clamp(0, width - 1);
    final z0 = gz.floor().clamp(0, height - 1);
    final x1 = math.min(width - 1, x0 + 1);
    final z1 = math.min(height - 1, z0 + 1);
    final tx = gx - x0;
    final tz = gz - z0;
    final h00 = surfaceHeightAt(x0, z0);
    final h10 = surfaceHeightAt(x1, z0);
    final h01 = surfaceHeightAt(x0, z1);
    final h11 = surfaceHeightAt(x1, z1);
    final a = h00 + (h10 - h00) * tx;
    final b = h01 + (h11 - h01) * tx;
    return a + (b - a) * tz;
  }

  /// Steps hydrostatic pipe fluxes, free-surface heights & velocities forward by [dt].
  void step(double dt) {
    if (dt <= 0.0 || width <= 0 || height <= 0) return;
    final nSub = math.max(1, (dt / 0.02).ceil());
    final h = dt / nSub;
    for (var s = 0; s < nSub; s++) {
      _stepSubstep(h);
    }
  }

  void _stepSubstep(double dt) {
    final area = math.max(1e-6, cellSize * cellSize);
    final pipeAccel = dt * pipeCrossSection * gravity / math.max(1e-4, cellSize);

    // 1. Update virtual outflow fluxes with mass-conserving capacity clamp
    for (var z = 0; z < height; z++) {
      for (var x = 0; x < width; x++) {
        final i = _idx(x, z);
        final d = _depth[i];
        if (d <= 1e-9) {
          _fluxL[i] = 0.0;
          _fluxR[i] = 0.0;
          _fluxU[i] = 0.0;
          _fluxD[i] = 0.0;
          continue;
        }
        final selfH = _terrain[i] + d;

        var fL = x > 0
            ? math.max(0.0, _fluxL[i] * damping + pipeAccel * (selfH - (_terrain[_idx(x - 1, z)] + _depth[_idx(x - 1, z)])))
            : 0.0;
        var fR = x + 1 < width
            ? math.max(0.0, _fluxR[i] * damping + pipeAccel * (selfH - (_terrain[_idx(x + 1, z)] + _depth[_idx(x + 1, z)])))
            : 0.0;
        var fU = z > 0
            ? math.max(0.0, _fluxU[i] * damping + pipeAccel * (selfH - (_terrain[_idx(x, z - 1)] + _depth[_idx(x, z - 1)])))
            : 0.0;
        var fD = z + 1 < height
            ? math.max(0.0, _fluxD[i] * damping + pipeAccel * (selfH - (_terrain[_idx(x, z + 1)] + _depth[_idx(x, z + 1)])))
            : 0.0;

        final outSum = fL + fR + fU + fD;
        if (outSum > 1e-12) {
          final k = math.min(1.0, (d * area) / (outSum * dt));
          fL *= k;
          fR *= k;
          fU *= k;
          fD *= k;
        }
        _fluxL[i] = fL;
        _fluxR[i] = fR;
        _fluxU[i] = fU;
        _fluxD[i] = fD;
      }
    }

    // 2. Apply net inflow - outflow volume exchange and derive velocity vectors
    for (var z = 0; z < height; z++) {
      for (var x = 0; x < width; x++) {
        final i = _idx(x, z);
        final inFromL = x > 0 ? _fluxR[_idx(x - 1, z)] : 0.0;
        final inFromR = x + 1 < width ? _fluxL[_idx(x + 1, z)] : 0.0;
        final inFromU = z > 0 ? _fluxD[_idx(x, z - 1)] : 0.0;
        final inFromD = z + 1 < height ? _fluxU[_idx(x, z + 1)] : 0.0;

        final outL = _fluxL[i];
        final outR = _fluxR[i];
        final outU = _fluxU[i];
        final outD = _fluxD[i];

        final netVolume = dt * ((inFromL + inFromR + inFromU + inFromD) - (outL + outR + outU + outD));
        final oldD = _depth[i];
        var newD = math.max(0.0, oldD + netVolume / area);
        if (evaporationRate > 0.0 && newD > 0.0) {
          newD = math.max(0.0, newD - evaporationRate * dt);
        }
        _nextDepth[i] = newD;

        final avgD = 0.5 * (oldD + newD);
        if (avgD > 1e-4) {
          final netQx = 0.5 * ((inFromL - outL) + (outR - inFromR));
          final netQz = 0.5 * ((inFromU - outU) + (outD - inFromD));
          _vx[i] = netQx / (cellSize * avgD);
          _vz[i] = netQz / (cellSize * avgD);
        } else {
          _vx[i] = 0.0;
          _vz[i] = 0.0;
        }
      }
    }

    for (var i = 0; i < _depth.length; i++) {
      _depth[i] = _nextDepth[i];
    }
  }
}
