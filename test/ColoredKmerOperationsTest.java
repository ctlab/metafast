import algo.ColoredKmerOperations;
import org.junit.Test;
import static org.junit.Assert.assertEquals;

public class ColoredKmerOperationsTest {

    private final int POWER = 20;
    private final long MAX_VALUE = (1L << POWER) - 1;

    private void show(long v) {
        System.out.println(Long.toBinaryString(v));
    }

    @Test
    public void testGetValue() {
        long a = 0;
        long b = 239;
        long c = MAX_VALUE;

        long x = a << 0 * POWER;
        long y = b << 1 * POWER;
        long z = c << 2 * POWER;

        assertEquals(a, ColoredKmerOperations.getValue(x, 0));
        assertEquals(b, ColoredKmerOperations.getValue(y, 1));
        assertEquals(c, ColoredKmerOperations.getValue(z, 2));

        long v = x + y + z;
        assertEquals(a, ColoredKmerOperations.getValue(v, 0));
        assertEquals(b, ColoredKmerOperations.getValue(v, 1));
        assertEquals(c, ColoredKmerOperations.getValue(v, 2));
    }

    @Test
    public void testAddOne() {
        long x = 0;

        x = ColoredKmerOperations.addValue(x, 0);
        x = ColoredKmerOperations.addValue(x, 1);
        x = ColoredKmerOperations.addValue(x, 2);

        assertEquals(1, ColoredKmerOperations.getValue(x, 0));
        assertEquals(1, ColoredKmerOperations.getValue(x, 1));
        assertEquals(1, ColoredKmerOperations.getValue(x, 2));

        x = ColoredKmerOperations.addValue(x, 0);
        x = ColoredKmerOperations.addValue(x, 1);
        x = ColoredKmerOperations.addValue(x, 2);
        x = ColoredKmerOperations.addValue(x, 0);
        x = ColoredKmerOperations.addValue(x, 1);
        x = ColoredKmerOperations.addValue(x, 2);

        assertEquals(3, ColoredKmerOperations.getValue(x, 0));
        assertEquals(3, ColoredKmerOperations.getValue(x, 1));
        assertEquals(3, ColoredKmerOperations.getValue(x, 2));

        long a = MAX_VALUE;

        long t = a << 0 * POWER;
        long y = a << 1 * POWER;
        long z = a << 2 * POWER;

        long v = t + y + z;
        v = ColoredKmerOperations.addValue(v, 0);
        v = ColoredKmerOperations.addValue(v, 1);
        v = ColoredKmerOperations.addValue(v, 2);

        assertEquals(MAX_VALUE, ColoredKmerOperations.getValue(v, 0));
        assertEquals(MAX_VALUE, ColoredKmerOperations.getValue(v, 1));
        assertEquals(MAX_VALUE, ColoredKmerOperations.getValue(v, 2));
    }

    @Test
    public void testAddValue() {
        long v = 0;

        short a = 1;
        short b = 239;
        short c = 0;

        v = ColoredKmerOperations.addValue(v, 0, a);
        v = ColoredKmerOperations.addValue(v, 1, b);
        v = ColoredKmerOperations.addValue(v, 2, c);

        assertEquals(a, ColoredKmerOperations.getValue(v, 0));
        assertEquals(b, ColoredKmerOperations.getValue(v, 1));
        assertEquals(c, ColoredKmerOperations.getValue(v, 2));


        long x = MAX_VALUE << 0 * POWER;
        long y = MAX_VALUE << 1 * POWER;
        long z = MAX_VALUE << 2 * POWER;

        v = x + y + z;

        assertEquals(MAX_VALUE, ColoredKmerOperations.getValue(v, 0));
        assertEquals(MAX_VALUE, ColoredKmerOperations.getValue(v, 1));
        assertEquals(MAX_VALUE, ColoredKmerOperations.getValue(v, 2));

        v = ColoredKmerOperations.addValue(v, 0, (short)0);
        v = ColoredKmerOperations.addValue(v, 1, (short)1);
        v = ColoredKmerOperations.addValue(v, 2, (short)239);

        assertEquals(MAX_VALUE, ColoredKmerOperations.getValue(v, 0));
        assertEquals(MAX_VALUE, ColoredKmerOperations.getValue(v, 1));
        assertEquals(MAX_VALUE, ColoredKmerOperations.getValue(v, 2));
    }

    @Test
    public void testSetValue() {
        long v = 0;

        long a = 1;
        long b = 239;
        long c = MAX_VALUE;

        v = ColoredKmerOperations.setValue(v, 0, a);
        v = ColoredKmerOperations.setValue(v, 1, b);
        v = ColoredKmerOperations.setValue(v, 2, c);

        assertEquals(a, ColoredKmerOperations.getValue(v, 0));
        assertEquals(b, ColoredKmerOperations.getValue(v, 1));
        assertEquals(c, ColoredKmerOperations.getValue(v, 2));

        a = 357;
        b = 0;
        c = 98765;

        v = ColoredKmerOperations.setValue(v, 0, a);
        v = ColoredKmerOperations.setValue(v, 1, b);
        v = ColoredKmerOperations.setValue(v, 2, c);

        assertEquals(a, ColoredKmerOperations.getValue(v, 0));
        assertEquals(b, ColoredKmerOperations.getValue(v, 1));
        assertEquals(c, ColoredKmerOperations.getValue(v, 2));

    }

    @Test
    public void testClearValue() {
        long a = 1;
        long b = 239;
        long c = MAX_VALUE;

        long x = a << 0 * POWER;
        long y = b << 1 * POWER;
        long z = c << 2 * POWER;

        long v = x + y + z;
        assertEquals(a, ColoredKmerOperations.getValue(v, 0));
        assertEquals(b, ColoredKmerOperations.getValue(v, 1));
        assertEquals(c, ColoredKmerOperations.getValue(v, 2));

        v = ColoredKmerOperations.clearValue(v, 0);
        assertEquals(0, ColoredKmerOperations.getValue(v, 0));
        assertEquals(b, ColoredKmerOperations.getValue(v, 1));
        assertEquals(c, ColoredKmerOperations.getValue(v, 2));

        v = ColoredKmerOperations.clearValue(v, 1);
        assertEquals(0, ColoredKmerOperations.getValue(v, 0));
        assertEquals(0, ColoredKmerOperations.getValue(v, 1));
        assertEquals(c, ColoredKmerOperations.getValue(v, 2));

        v = ColoredKmerOperations.clearValue(v, 2);
        assertEquals(0, ColoredKmerOperations.getValue(v, 0));
        assertEquals(0, ColoredKmerOperations.getValue(v, 1));
        assertEquals(0, ColoredKmerOperations.getValue(v, 2));

    }
}