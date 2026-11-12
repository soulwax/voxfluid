import 'package:test/test.dart';
import 'package:voxfluid/voxfluid.dart';

void main() {
  test('ShallowWaterGrid conserves fluid volume and spills over a weir ridge', () {
    final water = ShallowWaterGrid(12, 8, cellSize: 1.0, damping: 0.99);
    // Raise a weir wall at x == 5 with height 0.6m
    for (var z = 0; z < 8; z++) {
      water.setTerrainHeight(5, z, 0.6);
    }
    // Pour deep water column in upper basin (3*4*2.8 = 33.6 m^3 > 40*0.6 = 24 m^3 dam capacity)
    for (var z = 2; z <= 5; z++) {
      for (var x = 1; x <= 3; x++) {
        water.addWaterDepth(x, z, 2.8);
      }
    }
    final initialVolume = water.totalWaterVolume;
    expect(water.waterDepthAt(8, 4), closeTo(0.0, 1e-9));

    water.step(3.0);

    // Volume conserved without evaporation
    expect(water.totalWaterVolume, closeTo(initialVolume, 1e-4));
    // Crest overflows x == 5 ridge into lower basin (x >= 6)
    expect(water.waterDepthAt(8, 4), greaterThan(0.05));
    expect(water.floodedCellCount, greaterThan(20));
  });

  test('CellularFireGrid spreads downwind, generates smoke & extinguishes with water radius', () {
    final fire = CellularFireGrid(
      17,
      9,
      windX: 2.2,
      windZ: 0.0,
      spreadRate: 2.0,
    )..fillFuelUniform(1.0);

    fire.ignite(8, 4, heat: 2.4);
    final initialFuel = fire.totalFuelRemaining;

    fire.step(2.5);

    expect(fire.burningCellCount, greaterThan(3));
    expect(fire.totalSmoke, greaterThan(0.5));
    expect(fire.totalFuelRemaining, lessThan(initialFuel));

    // Downwind side (+x from center 8) should receive strictly more heat/burn than upwind (-x)
    var upwindHeat = 0.0;
    var downwindHeat = 0.0;
    for (var d = 1; d <= 5; d++) {
      upwindHeat += fire.heatAt(8 - d, 4);
      downwindHeat += fire.heatAt(8 + d, 4);
    }
    expect(downwindHeat, greaterThan(upwindHeat));

    // Quench whole area
    fire.extinguishRadius(8, 4, 12.0, waterStrength: 2.0);
    expect(fire.burningCellCount, 0);
  });

  test('ShallowWaterGrid flooding extinguishes coupled CellularFireGrid cells', () {
    final water = ShallowWaterGrid(10, 6);
    final fire = CellularFireGrid(10, 6)..fillFuelUniform(1.0);
    fire.ignite(5, 3, heat: 2.2);
    expect(fire.burningCellCount, 1);

    water.addWaterVolume(2, 3, 18.0);
    for (var s = 0; s < 30; s++) {
      water.step(0.1);
      fire.step(0.1);
      for (var z = 0; z < 6; z++) {
        for (var x = 0; x < 10; x++) {
          if (water.waterDepthAt(x, z) > 0.03) {
            fire.quenchCell(x, z, moistureAdd: 1.0);
          }
        }
      }
    }

    expect(water.waterDepthAt(5, 3), greaterThan(0.03));
    expect(fire.burnAt(5, 3), closeTo(0.0, 1e-6));
  });
}
