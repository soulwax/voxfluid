import 'dart:math' as math;
import 'dart:typed_data';

/// Cellular 2D combustion, wind-advected fire spread & smoke dissipation grid.
final class CellularFireGrid {
  final int width;
  final int height;
  final double cellSize;
  double windX;
  double windZ;
  double ignitionThreshold;
  double coolingRate;
  double spreadRate;
  double burnRate;
  double smokeProductionRate;
  double smokeDecayRate;

  final Float64List _fuel;
  final Float64List _heat;
  final Float64List _nextHeat;
  final Float64List _burn;
  final Float64List _smoke;
  final Float64List _moisture;

  CellularFireGrid(
    this.width,
    this.height, {
    this.cellSize = 1.0,
    this.windX = 0.0,
    this.windZ = 0.0,
    this.ignitionThreshold = 0.75,
    this.coolingRate = 0.35,
    this.spreadRate = 1.8,
    this.burnRate = 0.22,
    this.smokeProductionRate = 1.3,
    this.smokeDecayRate = 0.45,
  })  : _fuel = Float64List(width * height),
        _heat = Float64List(width * height),
        _nextHeat = Float64List(width * height),
        _burn = Float64List(width * height),
        _smoke = Float64List(width * height),
        _moisture = Float64List(width * height);

  int _idx(int x, int z) => z * width + x;
  bool inBounds(int x, int z) => x >= 0 && x < width && z >= 0 && z < height;

  void setCell(int x, int z, {double fuel = 1.0, double moisture = 0.0}) {
    if (!inBounds(x, z)) return;
    final i = _idx(x, z);
    _fuel[i] = math.max(0.0, fuel);
    _moisture[i] = math.max(0.0, moisture);
  }

  void fillFuelUniform(double fuel, {double moisture = 0.0}) {
    _fuel.fillRange(0, _fuel.length, math.max(0.0, fuel));
    _moisture.fillRange(0, _moisture.length, math.max(0.0, moisture));
  }

  void ignite(int x, int z, {double heat = 2.0}) {
    if (!inBounds(x, z)) return;
    final i = _idx(x, z);
    if (_fuel[i] <= 1e-4) return;
    _heat[i] = math.max(_heat[i], heat);
    if (_moisture[i] <= 0.1) {
      _burn[i] = math.max(_burn[i], 0.8);
    }
  }

  void quenchCell(int x, int z, {double moistureAdd = 0.6}) {
    if (!inBounds(x, z)) return;
    final i = _idx(x, z);
    _burn[i] = 0.0;
    _heat[i] = 0.0;
    _moisture[i] = math.min(2.0, _moisture[i] + moistureAdd);
  }

  void extinguishRadius(int cx, int cz, double radiusCells, {double waterStrength = 1.0}) {
    final rCeil = radiusCells.ceil();
    final r2 = radiusCells * radiusCells;
    for (var dz = -rCeil; dz <= rCeil; dz++) {
      for (var dx = -rCeil; dx <= rCeil; dx++) {
        if (dx * dx + dz * dz > r2) continue;
        quenchCell(cx + dx, cz + dz, moistureAdd: waterStrength * 0.6);
      }
    }
  }

  double fuelAt(int x, int z) => inBounds(x, z) ? _fuel[_idx(x, z)] : 0.0;
  double heatAt(int x, int z) => inBounds(x, z) ? _heat[_idx(x, z)] : 0.0;
  double burnAt(int x, int z) => inBounds(x, z) ? _burn[_idx(x, z)] : 0.0;
  double smokeAt(int x, int z) => inBounds(x, z) ? _smoke[_idx(x, z)] : 0.0;
  double moistureAt(int x, int z) => inBounds(x, z) ? _moisture[_idx(x, z)] : 0.0;

  int get burningCellCount {
    var count = 0;
    for (var i = 0; i < _burn.length; i++) {
      if (_burn[i] > 0.05) count++;
    }
    return count;
  }

  double get totalFuelRemaining {
    var sum = 0.0;
    for (var i = 0; i < _fuel.length; i++) {
      sum += _fuel[i];
    }
    return sum;
  }

  double get totalSmoke {
    var sum = 0.0;
    for (var i = 0; i < _smoke.length; i++) {
      sum += _smoke[i];
    }
    return sum;
  }

  void step(double dt) {
    if (dt <= 0.0 || width <= 0 || height <= 0) return;
    final nSub = math.max(1, (dt / 0.05).ceil());
    final h = dt / nSub;
    for (var s = 0; s < nSub; s++) {
      _stepSubstep(h);
    }
  }

  void _stepSubstep(double dt) {
    for (var i = 0; i < _heat.length; i++) {
      _nextHeat[i] = _heat[i];
    }

    const dirs = [
      (-1, 0, 1.0),
      (1, 0, 1.0),
      (0, -1, 1.0),
      (0, 1, 1.0),
      (-1, -1, 0.707),
      (1, -1, 0.707),
      (-1, 1, 0.707),
      (1, 1, 0.707),
    ];

    for (var z = 0; z < height; z++) {
      for (var x = 0; x < width; x++) {
        final i = _idx(x, z);
        final burnVal = _burn[i];
        if (burnVal <= 0.02) continue;

        final emitted = burnVal * spreadRate * dt;
        _nextHeat[i] = math.min(4.0, _nextHeat[i] + burnVal * 1.1 * dt);
        for (final (dx, dz, baseW) in dirs) {
          final nx = x + dx;
          final nz = z + dz;
          if (!inBounds(nx, nz)) continue;
          final ni = _idx(nx, nz);
          if (_fuel[ni] <= 1e-4) continue;
          final len = math.sqrt(dx * dx + dz * dz);
          final windAlign = (dx * windX + dz * windZ) / len;
          final windBias = (1.0 + windAlign * 0.45).clamp(0.15, 3.0);
          _nextHeat[ni] = math.min(4.0, _nextHeat[ni] + emitted * baseW * windBias);
        }
      }
    }

    for (var i = 0; i < _heat.length; i++) {
      var h = math.max(0.0, _nextHeat[i] - coolingRate * dt);
      if (_moisture[i] > 0.0) {
        final dry = math.min(_moisture[i], h * 0.65 * dt);
        _moisture[i] -= dry;
        h = math.max(0.0, h - dry * 1.4);
      }

      if (_fuel[i] > 1e-4) {
        if (h >= ignitionThreshold && _moisture[i] <= 0.05) {
          _burn[i] = math.min(1.5, _burn[i] + (h - ignitionThreshold + 0.5) * 2.2 * dt);
        } else {
          _burn[i] = math.max(0.0, _burn[i] - 0.9 * dt);
        }
        if (_burn[i] > 0.02) {
          final burned = math.min(_fuel[i], _burn[i] * burnRate * dt);
          _fuel[i] -= burned;
          _smoke[i] = math.min(3.0, _smoke[i] + _burn[i] * smokeProductionRate * dt);
          if (_fuel[i] <= 1e-4) {
            _fuel[i] = 0.0;
            _burn[i] = 0.0;
          }
        }
      } else {
        _burn[i] = 0.0;
      }

      _heat[i] = h;
      _smoke[i] = math.max(0.0, _smoke[i] - smokeDecayRate * dt);
    }
  }
}
