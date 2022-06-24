package algo;

public class ColoredKmerOperations {

    private static final int POWER = 20;

    public static long getValue(long value, int color) {
        long mask = (1L<<(color+1)*POWER) - (1L<<color*POWER);
        return (value&mask) >> (color*POWER);
    }

    public static long addValue(long value, int color) {
        return addValue(value, color, (short)1);
    }

    public static long addValue(long value, int color, short add) {
        long colorValue = getValue(value, color);
        long newColorValue = Math.min(colorValue + add, (1L<<POWER) - 1);
        long newValue = setValue(value, color, newColorValue);
        return newValue;
    }

    public static long setValue(long value, int color, long newColorValue) {
        long newValue = clearValue(value, color);
        newValue = newValue | (newColorValue << (color*POWER));
        return newValue;
    }

    public static long clearValue(long value, int color) {
        long mask = (1L<<(color+1)*POWER) - (1L<<color*POWER);
        mask = ~mask;
        long newValue = value&mask;
        return newValue;
    }

}
