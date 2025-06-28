package fragrant;

public class EndUtil {

    public static int getHeight(long seed, int x, int z) {
        return BedrockEndTerrainGen.getEndSurfaceHeight(seed, x, z);
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