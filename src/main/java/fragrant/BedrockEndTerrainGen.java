package fragrant;

public class BedrockEndTerrainGen {

    public static class EndGenerator {
        private final EndNoise en;
        private final SurfaceNoise sn;

        public EndGenerator(long seed) {
            BedrockRandom mt = new BedrockRandom(seed);
            mt.skipNextN(17292);
            en = new EndNoise(mt);

            mt = new BedrockRandom(seed);
            sn = new SurfaceNoise(mt);
        }

        public double getDensity(int x, int y, int z) {
            int cx = x >> 3, cz = z >> 3, cy = y >> 2;
            double dx = (x & 7) * 0.125, dz = (z & 7) * 0.125, dy = (y & 3) * 0.25;

            double[] c = new double[8];
            sampleCol(c, 0, cx, cz, cy);
            sampleCol(c, 2, cx, cz + 1, cy);
            sampleCol(c, 4, cx + 1, cz, cy);
            sampleCol(c, 6, cx + 1, cz + 1, cy);

            return lerp3D(dy, dx, dz, c[0], c[1], c[4], c[5], c[2], c[3], c[6], c[7]);
        }

        private void sampleCol(double[] c, int idx, int x, int z, int y) {
            double[] up = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,63./64,62./64,61./64,60./64,59./64,58./64,57./64,56./64,55./64,54./64,53./64,52./64,51./64,50./64,49./64,48./64,47./64,46./64};
            double[] lo = {0,0,1./7,2./7,3./7,4./7,5./7,6./7,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
            double d = endHeight(x, z) - 8;
            for (int i = 0; i < 2; i++) {
                int yy = y + i;
                if (lo[yy] == 0) {
                    c[idx + i] = -30;
                } else {
                    double n = surfaceNoise(x, yy, z, -128, 128) + d;
                    n = lerp(up[yy], -3000, n);
                    c[idx + i] = lerp(lo[yy], -30, n);
                }
            }
        }

        private double endHeight(int x, int z) {
            int hx = x / 2, hz = z / 2, ox = x % 2, oz = z % 2;
            long h = 64L * (x * (long) x + z * (long) z);

            for (int j = -12; j <= 12; j++) {
                for (int i = -12; i <= 12; i++) {
                    long rx = hx + i, rz = hz + j;
                    if (rx * rx + rz * rz > 4096 && simplex2D(rx, rz) < -0.9) {
                        int v = (int) ((Math.abs(rx) * 3439 + Math.abs(rz) * 147) % 13 + 9);
                        rx = ox - i * 2;
                        rz = oz - j * 2;
                        long n = (rx * rx + rz * rz) * v * v;
                        if (n < h) h = n;
                    }
                }
            }
            return Math.max(-100, Math.min(80, 100 - (float) Math.sqrt(h)));
        }

        private double simplex2D(double x, double y) {
            final double S = 0.366025403784;
            final double U = 0.211324865405;

            double h = (x + y) * S;
            int hx = (int) Math.floor(x + h);
            int hz = (int) Math.floor(y + h);

            double m = (hx + hz) * U;
            double x0 = x - (hx - m), y0 = y - (hz - m);

            int ox = x0 > y0 ? 1 : 0, oz = 1 - ox;
            double x1 = x0 - ox + U, y1 = y0 - oz + U;
            double x2 = x0 - 1 + 2 * U, y2 = y0 - 1 + 2 * U;

            int[] d = en.d;
            int g0 = d[0xff & (d[0xff & hz] + hx)];
            int g1 = d[0xff & (d[0xff & (hz + oz)] + hx + ox)];
            int g2 = d[0xff & (d[0xff & (hz + 1)] + hx + 1)];

            return 70 * (simplexGrad(g0 % 12, x0, y0) + simplexGrad(g1 % 12, x1, y1) + simplexGrad(g2 % 12, x2, y2));
        }

        private double simplexGrad(int i, double x, double y) {
            double c = 0.5 - x * x - y * y;
            return c < 0 ? 0 : c * c * c * c * grad(i, x, y, 0);
        }

        private double surfaceNoise(int x, int y, int z, double nmin, double nmax) {
            double xzs = 1368.824, ys = 684.412;
            double vmin = 0, vmax = 0, p = 1.0 / 32768, a = 64;

            for (int i = 15; i >= 0; i--, a *= 0.5, p *= 2) {
                double dx = x * xzs * p, dz = z * xzs * p, dy = y * ys * p;
                vmin += perlin(sn.min[i], dx, dy, dz, ys * p, dy) * a;
                vmax += perlin(sn.max[i], dx, dy, dz, ys * p, dy) * a;
                if (vmin - a > nmax && vmax - a > nmax) return nmax;
                if (vmin + a < nmin && vmax + a < nmin) return nmin;
            }

            double xs = 17.1103, yst = 4.27758, vm = 0.5;
            p = 1.0 / 128;
            a = 6.4;

            for (int i = 7; i >= 0; i--, a *= 0.5, p *= 2) {
                double dx = x * xs * p, dz = z * xs * p, dy = y * yst * p;
                vm += perlin(sn.main[i], dx, dy, dz, yst * p, dy) * a;
                if (vm - a > 1) return vmax;
                if (vm + a < 0) return vmin;
            }

            return vm <= 0 ? vmin : vm >= 1 ? vmax : vmin + vm * (vmax - vmin);
        }

        private double perlin(Perlin p, double d1, double d2, double d3, double ya, double ym) {
            int h1, h2, h3;
            double t1, t2, t3;

            if (d2 == 0) {
                d2 = p.d2;
                h2 = p.h2;
                t2 = p.t2;
            } else {
                d2 += p.b;
                double i2 = Math.floor(d2);
                d2 -= i2;
                h2 = (int) i2;
                t2 = d2 * d2 * d2 * (d2 * (d2 * 6 - 15) + 10);
            }

            d1 += p.a;
            d3 += p.c;
            double i1 = Math.floor(d1), i3 = Math.floor(d3);
            d1 -= i1;
            d3 -= i3;
            h1 = (int) i1;
            h3 = (int) i3;
            t1 = d1 * d1 * d1 * (d1 * (d1 * 6 - 15) + 10);
            t3 = d3 * d3 * d3 * (d3 * (d3 * 6 - 15) + 10);

            if (ya != 0) {
                d2 -= Math.floor(Math.min(ym, d2) / ya) * ya;
            }

            int[] d = p.d;
            int a1 = d[h1 & 0xff] + h2, b1 = d[(h1 + 1) & 0xff] + h2;
            int a2 = d[a1 & 0xff] + h3, b2 = d[b1 & 0xff] + h3;
            int a3 = d[(a1 + 1) & 0xff] + h3, b3 = d[(b1 + 1) & 0xff] + h3;

            double l1 = grad(d[a2 & 0xff], d1, d2, d3);
            double l2 = grad(d[b2 & 0xff], d1 - 1, d2, d3);
            double l3 = grad(d[a3 & 0xff], d1, d2 - 1, d3);
            double l4 = grad(d[b3 & 0xff], d1 - 1, d2 - 1, d3);
            double l5 = grad(d[(a2 + 1) & 0xff], d1, d2, d3 - 1);
            double l6 = grad(d[(b2 + 1) & 0xff], d1 - 1, d2, d3 - 1);
            double l7 = grad(d[(a3 + 1) & 0xff], d1, d2 - 1, d3 - 1);
            double l8 = grad(d[(b3 + 1) & 0xff], d1 - 1, d2 - 1, d3 - 1);

            return lerp(t3, lerp(t2, lerp(t1, l1, l2), lerp(t1, l3, l4)), lerp(t2, lerp(t1, l5, l6), lerp(t1, l7, l8)));
        }
    }

    static class Perlin {
        double a, b, c, d2, t2;
        int[] d = new int[257];
        int h2;

        Perlin(BedrockRandom rng) {
            a = rng.nextDouble() * 256;
            b = rng.nextDouble() * 256;
            c = rng.nextDouble() * 256;

            for (int i = 0; i < 256; i++) d[i] = i;
            for (int i = 0; i < 256; i++) {
                int j = rng.nextInt(256 - i) + i;
                int t = d[i];
                d[i] = d[j];
                d[j] = t;
            }
            d[256] = d[0];

            double i2 = Math.floor(b);
            d2 = b - i2;
            h2 = (int) i2;
            t2 = d2 * d2 * d2 * (d2 * (d2 * 6 - 15) + 10);
        }
    }

    static class EndNoise {
        int[] d;
        EndNoise(BedrockRandom rng) {
            Perlin p = new Perlin(rng);
            d = p.d;
        }
    }

    static class SurfaceNoise {
        Perlin[] min = new Perlin[16];
        Perlin[] max = new Perlin[16];
        Perlin[] main = new Perlin[8];

        SurfaceNoise(BedrockRandom rng) {
            for (int i = 0; i < 40; i++) {
                Perlin p = new Perlin(rng);
                if (i < 16) min[i] = p;
                else if (i < 32) max[i - 16] = p;
                else main[i - 32] = p;
            }
        }
    }

    private static double grad(int i, double x, double y, double z) {
        return switch (i & 15) {
            case 0 -> x + y;
            case 1 -> -x + y;
            case 2 -> x - y;
            case 3 -> -x - y;
            case 4 -> x + z;
            case 5 -> -x + z;
            case 6 -> x - z;
            case 7 -> -x - z;
            case 8 -> y + z;
            case 9 -> -y + z;
            case 10 -> y - z;
            case 11 -> -y - z;
            case 12 -> x + y;
            case 13 -> -y + z;
            case 14 -> -x + y;
            case 15 -> -y - z;
            default -> 0;
        };
    }

    private static double lerp(double t, double a, double b) {
        return a + t * (b - a);
    }

    private static double lerp3D(double dx, double dy, double dz, double v000, double v100, double v010, double v110, double v001, double v101, double v011, double v111) {
        return lerp(dz, lerp(dy, lerp(dx, v000, v100), lerp(dx, v010, v110)), lerp(dy, lerp(dx, v001, v101), lerp(dx, v011, v111)));
    }

    public static int getEndSurfaceHeight(long seed, int x, int z) {
        EndGenerator gen = new EndGenerator(seed);
        for (int y = 127; y >= 0; y--) {
            if (gen.getDensity(x, y, z) > 0) return y;
        }
        return 0;
    }
}