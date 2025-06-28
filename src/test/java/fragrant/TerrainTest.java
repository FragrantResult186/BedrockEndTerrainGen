package fragrant;

import static org.junit.jupiter.api.Assertions.*;

public class TerrainTest {

    @org.junit.jupiter.api.Test
    public void testHeight() {
        assertEquals(0, EndUtil.getHeight(93376L, 0, 0));
        assertEquals(59, EndUtil.getHeight(-99999999999L, 0, 0));
    }

    @org.junit.jupiter.api.Test
    public void testBlock() {
        assertTrue(EndUtil.isEndStone(560233420133338179L, -130, 55, -2));
        assertTrue(EndUtil.isAir(560233420133338179L, -130, 54, -2));
    }

}
