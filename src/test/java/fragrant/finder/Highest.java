package fragrant.finder;

import fragrant.EndUtil;

public class Highest {
    public static void main(String[] args) {
        int min = 0;
        for (long seed = 13541L; seed < Integer.MAX_VALUE; seed++) {
            int y = EndUtil.getHeight(seed, 0, 0);
            if (y >= min) {
                min = y;
                System.out.printf("%d %d\n", seed, y);
            }
        }
    }
}
