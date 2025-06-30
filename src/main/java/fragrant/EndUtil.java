package fragrant;

public class EndUtil {

    public static int getHeight(long seed, int x, int z) {
        return BedrockEndTerrainGen.getEndSurfaceHeight(seed, x, z);
    }

    public static int[][] getHeightMap(long seed, int minX, int minZ, int maxX, int maxZ) {
        int width = maxX - minX + 1;
        int height = maxZ - minZ + 1;
        int[][] heightMap = new int[height][width];
        java.util.stream.IntStream.range(0, height).parallel().forEach(zIndex -> {
            int z = minZ + zIndex;
            BedrockEndTerrainGen.EndGenerator gen = new BedrockEndTerrainGen.EndGenerator(seed);
            for (int x = minX; x <= maxX; x++) {
                for (int y = 80; y >= 0; y--) {
                    if (gen.getDensity(x, y, z) > 0) {
                        heightMap[zIndex][x - minX] = y;
                        break;
                    }
                }
            }
        });
        return heightMap;
    }

    public static double getDensity(long seed, int x, int y, int z) {
        BedrockEndTerrainGen.EndGenerator gen = new BedrockEndTerrainGen.EndGenerator(seed);
        return gen.getDensity(x, y, z);
    }

    public static boolean isAir(long seed, int x, int y, int z) {
        return getDensity(seed, x, y, z) <= 0;
    }

    public static boolean isEndStone(long seed, int x, int y, int z) {
        return getDensity(seed, x, y, z) > 0;
    }

}