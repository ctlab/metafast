package structures.map;

import it.unimi.dsi.fastutil.HashCommon;
import org.junit.Before;
import org.junit.Test;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import static org.junit.Assert.*;

public class Long2BitLongaHashMapTest {
    Random rand;

    @Before
    public void before() {
        rand = new Random(239);
    }
    
    @Test
    public void testCreate() {
        Long2BitLongaHashMap hm = new Long2BitLongaHashMap(20, 0.75f, 130);
        assertEquals(130, hm.data.sizeBitSet);
        assertEquals(3, hm.data.valueForFreeKey.length);
        assertNull(hm.data.values[0]);

        hm.reset();

        assertEquals(130, hm.data.sizeBitSet);
        assertEquals(3, hm.data.valueForFreeKey.length);
        assertNull(hm.data.values[0]);
    }

    @Test
    public void testAdd() {
        Long2BitLongaHashMap hm = new Long2BitLongaHashMap(20, 0.75f, 130);
        assertNull(hm.data.values[getPositionInt(hm.data, 12)]);

        hm.add(12);
        assertEquals(1, hm.data.values[getPositionInt(hm.data, 12)][0]);

        hm.add(12);
        assertEquals(1, hm.data.values[getPositionInt(hm.data, 12)][0]);

        hm.reset();
        assertNull(hm.data.values[getPositionInt(hm.data, 12)]);

        int n_tests = 100000;
        int[] vals = new int[n_tests];
                for (int i = 0; i < n_tests; i++) {
            vals[i] = rand.nextInt(10000) + 1;
            hm.add(vals[i]);
        }

        for (int i = 0; i < n_tests; i++) {
            assertEquals(1, hm.data.values[getPositionInt(hm.data, vals[i])][0]);
        }

        hm.reset();

        for (int i = 0; i < n_tests; i++) {
            assertNull(hm.data.values[getPositionInt(hm.data, vals[i])]);
        }
    }

    @Test
    public void testSet() {
        Long2BitLongaHashMap hm = new Long2BitLongaHashMap(20, 0.75f, 130);

        long[] prev = hm.set(12, 62);
        long[] cur = new long[3];
        cur[0] = 1L<<62;
        assertNull(prev);
        assertArrayEquals(cur, hm.data.values[getPositionInt(hm.data, 12)]);

        prev = hm.set(12, 63);
        assertArrayEquals(cur, prev);
        cur[0] |= 1L<<63;
        assertArrayEquals(cur, hm.data.values[getPositionInt(hm.data, 12)]);

        prev = hm.set(12, 64);
        assertArrayEquals(cur, prev);
        cur[1] |= 1L<<0;
        assertArrayEquals(cur, hm.data.values[getPositionInt(hm.data, 12)]);

        prev = hm.set(12, 65);
        assertArrayEquals(cur, prev);
        cur[1] |= 1L<<1;
        assertArrayEquals(cur, hm.data.values[getPositionInt(hm.data, 12)]);

        hm.reset();
        assertNull(hm.data.values[getPositionInt(hm.data, 12)]);


        int n_tests = 10000;
        long[][] vals = new long[n_tests][];
                int[] keys = new int[n_tests];
        for (int i = 0; i < n_tests; i++) {
            boolean f = true;
            while (f) {
                keys[i] = rand.nextInt(100000) + 1;
                f = false;
                for (int j = 0; j < i; j++) {
                    if (keys[i] == keys[j]) {
                        f = true;
                        break;
                    }
                }
            }
            int n_sets = rand.nextInt(1000) + 1;
            for (int j = 0; j < n_sets; j++) {
                int val = rand.nextInt(130);
                prev = hm.set(keys[i], val);
                if (vals[i] == null) {
                    vals[i] = new long[3];
                    assertNull(prev);
                } else {
                    assertArrayEquals(vals[i], prev);
                }
                vals[i][val/64] |= 1L<<(val%64);
                assertArrayEquals(vals[i], hm.data.values[getPositionInt(hm.data, keys[i])]);
            }
        }

        for (int i = 0; i < n_tests; i++) {
            assertArrayEquals(vals[i], hm.data.values[getPositionInt(hm.data, keys[i])]);
        }

        hm.reset();

        for (int i = 0; i < n_tests; i++) {
            assertNull(hm.data.values[getPositionInt(hm.data, keys[i])]);
        }
    }

    @Test
    public void testGet() {
        Long2BitLongaHashMap hm = new Long2BitLongaHashMap(20, 0.75f, 130);

        assertNull(hm.get(12));

        hm.set(12, 62);
        long[] cur = new long[3];
        cur[0] = 1L<<62;
        assertArrayEquals(cur, hm.get(12));

        hm.set(12, 63);
        cur[0] |= 1L<<63;
        assertArrayEquals(cur, hm.get(12));

        hm.set(12, 64);
        cur[1] |= 1L<<0;
        assertArrayEquals(cur, hm.get(12));

        hm.set(12, 65);
        cur[1] |= 1L<<1;
        assertArrayEquals(cur, hm.get(12));

        hm.reset();
        assertNull(hm.get(12));


        int n_tests = 10000;
        long[][] vals = new long[n_tests][];
                int[] keys = new int[n_tests];
        for (int i = 0; i < n_tests; i++) {
            boolean f = true;
            while (f) {
                keys[i] = rand.nextInt(100000) + 1;
                f = false;
                for (int j = 0; j < i; j++) {
                    if (keys[i] == keys[j]) {
                        f = true;
                        break;
                    }
                }
            }
            int n_sets = rand.nextInt(1000) + 1;
            assertNull(hm.get(keys[i]));
            for (int j = 0; j < n_sets; j++) {
                int val = rand.nextInt(130);
                hm.set(keys[i], val);
                if (vals[i] == null) {
                    vals[i] = new long[3];
                }
                vals[i][val/64] |= 1L<<(val%64);
                assertArrayEquals(vals[i], hm.get(keys[i]));
            }
        }

        for (int i = 0; i < n_tests; i++) {
            assertArrayEquals(vals[i], hm.get(keys[i]));
        }

        hm.reset();

        for (int i = 0; i < n_tests; i++) {
            assertNull(hm.get(keys[i]));
        }
    }

    @Test
    public void testGetBitIndex() {
        Long2BitLongaHashMap hm = new Long2BitLongaHashMap(20, 0.75f, 130);

        for (int i = 0; i < 130; i++) {
            assertFalse(hm.get(12, i));
        }

        hm.set(12, 62);
        Set<Integer> ids = new HashSet<>();
        ids.add(62);
        for (int i = 0; i < 130; i++) {
            if (ids.contains(i)) {
                assertTrue(hm.get(12, i));
            } else {
                assertFalse(hm.get(12, i));
            }
        }


        hm.set(12, 63);
        ids.add(63);
        for (int i = 0; i < 130; i++) {
            if (ids.contains(i)) {
                assertTrue(hm.get(12, i));
            } else {
                assertFalse(hm.get(12, i));
            }
        }

        hm.set(12, 64);
        ids.add(64);
        for (int i = 0; i < 130; i++) {
            if (ids.contains(i)) {
                assertTrue(hm.get(12, i));
            } else {
                assertFalse(hm.get(12, i));
            }
        }

        hm.set(12, 65);
        ids.add(65);
        for (int i = 0; i < 130; i++) {
            if (ids.contains(i)) {
                assertTrue(hm.get(12, i));
            } else {
                assertFalse(hm.get(12, i));
            }
        }

        hm.reset();
        for (int i = 0; i < 130; i++) {
            assertFalse(hm.get(12, i));
        }


        int n_tests = 10000;
        Set<Integer>[] vals = new Set[n_tests];
                int[] keys = new int[n_tests];
        for (int i = 0; i < n_tests; i++) {
            boolean f = true;
            while (f) {
                keys[i] = rand.nextInt(100000) + 1;
                f = false;
                for (int j = 0; j < i; j++) {
                    if (keys[i] == keys[j]) {
                        f = true;
                        break;
                    }
                }
            }
            int n_sets = rand.nextInt(1000) + 1;
            for (int k = 0; k < 130; k++) {
                assertFalse(hm.get(keys[i], k));
            }
            for (int j = 0; j < n_sets; j++) {
                int val = rand.nextInt(130);
                hm.set(keys[i], val);
                if (vals[i] == null) {
                    vals[i] = new HashSet<>();
                }
                vals[i].add(val);
                for (int k = 0; k < 130; k++) {
                    if (vals[i].contains(k)) {
                        assertTrue(hm.get(keys[i], k));
                    } else {
                        assertFalse(hm.get(keys[i], k));
                    }
                }
            }
        }

        for (int i = 0; i < n_tests; i++) {
            for (int k = 0; k < 130; k++) {
                if (vals[i].contains(k)) {
                    assertTrue(hm.get(keys[i], k));
                } else {
                    assertFalse(hm.get(keys[i], k));
                }
            }
        }

        hm.reset();

        for (int i = 0; i < n_tests; i++) {
            for (int k = 0; k < 130; k++) {
                assertFalse(hm.get(keys[i], k));
            }
        }

    }

    @Test
    public void testGetWithEmpty() {
        Long2BitLongaHashMap hm = new Long2BitLongaHashMap(20, 0.75f, 130);
        long[] empty = new long[1];

        assertArrayEquals(empty, hm.getWithEmpty(12));

        hm.set(12, 62);
        long[] cur = new long[3];
        cur[0] = 1L<<62;
        assertArrayEquals(cur, hm.getWithEmpty(12));

        hm.set(12, 63);
        cur[0] |= 1L<<63;
        assertArrayEquals(cur, hm.getWithEmpty(12));

        hm.set(12, 64);
        cur[1] |= 1L<<0;
        assertArrayEquals(cur, hm.getWithEmpty(12));

        hm.set(12, 65);
        cur[1] |= 1L<<1;
        assertArrayEquals(cur, hm.getWithEmpty(12));

        hm.reset();
        assertArrayEquals(empty, hm.getWithEmpty(12));


        int n_tests = 10000;
        long[][] vals = new long[n_tests][];
                int[] keys = new int[n_tests];
        for (int i = 0; i < n_tests; i++) {
            boolean f = true;
            while (f) {
                keys[i] = rand.nextInt(100000) + 1;
                f = false;
                for (int j = 0; j < i; j++) {
                    if (keys[i] == keys[j]) {
                        f = true;
                        break;
                    }
                }
            }
            int n_sets = rand.nextInt(1000) + 1;
            assertArrayEquals(empty, hm.getWithEmpty(keys[i]));
            for (int j = 0; j < n_sets; j++) {
                int val = rand.nextInt(130);
                hm.set(keys[i], val);
                if (vals[i] == null) {
                    vals[i] = new long[3];
                }
                vals[i][val/64] |= 1L<<(val%64);
                assertArrayEquals(vals[i], hm.getWithEmpty(keys[i]));
            }
        }

        for (int i = 0; i < n_tests; i++) {
            assertArrayEquals(vals[i], hm.getWithEmpty(keys[i]));
        }

        hm.reset();

        for (int i = 0; i < n_tests; i++) {
            assertArrayEquals(empty, hm.getWithEmpty(keys[i]));
        }
    }

    @Test
    public void testGetCardinality() {
        Long2BitLongaHashMap hm = new Long2BitLongaHashMap(20, 0.75f, 130);

        assertEquals(0, hm.getCardinality(12));

        hm.set(12, 62);
        Set<Integer> ids = new HashSet<>();
        ids.add(62);
        assertEquals(ids.size(), hm.getCardinality(12));

        hm.set(12, 63);
        ids.add(63);
        assertEquals(ids.size(), hm.getCardinality(12));

        hm.set(12, 64);
        ids.add(64);
        assertEquals(ids.size(), hm.getCardinality(12));

        hm.set(12, 65);
        ids.add(65);
        assertEquals(ids.size(), hm.getCardinality(12));

        hm.reset();
        assertEquals(0, hm.getCardinality(12));


        int n_tests = 10000;
        Set<Integer>[] vals = new Set[n_tests];
                int[] keys = new int[n_tests];
        for (int i = 0; i < n_tests; i++) {
            boolean f = true;
            while (f) {
                keys[i] = rand.nextInt(100000) + 1;
                f = false;
                for (int j = 0; j < i; j++) {
                    if (keys[i] == keys[j]) {
                        f = true;
                        break;
                    }
                }
            }
            int n_sets = rand.nextInt(1000) + 1;
            assertEquals(0, hm.getCardinality(keys[i]));
            for (int j = 0; j < n_sets; j++) {
                int val = rand.nextInt(130);
                hm.set(keys[i], val);
                if (vals[i] == null) {
                    vals[i] = new HashSet<>();
                }
                vals[i].add(val);
                assertEquals(vals[i].size(), hm.getCardinality(keys[i]));
            }
        }

        for (int i = 0; i < n_tests; i++) {
            assertEquals(vals[i].size(), hm.getCardinality(keys[i]));
        }

        hm.reset();

        for (int i = 0; i < n_tests; i++) {
            assertEquals(0, hm.getCardinality(keys[i]));
        }
    }

    @Test
    public void testGetCardinalityRange() {
        Long2BitLongaHashMap hm = new Long2BitLongaHashMap(20, 0.75f, 130);

        assertEquals(0, hm.getCardinality(12, 62, 66));

        hm.set(12, 62);
        assertEquals(1, hm.getCardinality(12, 62, 63));

        hm.set(12, 63);
        assertEquals(1, hm.getCardinality(12, 63, 64));
        assertEquals(2, hm.getCardinality(12, 62, 64));

        hm.set(12, 64);
        assertEquals(1, hm.getCardinality(12, 64, 65));
        assertEquals(3, hm.getCardinality(12, 62, 65));

        hm.set(12, 65);
        assertEquals(1, hm.getCardinality(12, 65, 66));
        assertEquals(4, hm.getCardinality(12, 62, 66));

        hm.reset();
        assertEquals(0, hm.getCardinality(12, 62, 66));


        int n_tests = 10000;
        Set<Integer>[] vals = new Set[n_tests];
                int[] keys = new int[n_tests];
        for (int i = 0; i < n_tests; i++) {
            boolean f = true;
            while (f) {
                keys[i] = rand.nextInt(100000) + 1;
                f = false;
                for (int j = 0; j < i; j++) {
                    if (keys[i] == keys[j]) {
                        f = true;
                        break;
                    }
                }
            }
            int n_sets = rand.nextInt(1000) + 1;
            assertEquals(0, hm.getCardinality(keys[i], 1, 100));
            for (int j = 0; j < n_sets; j++) {
                int val = rand.nextInt(130);
                hm.set(keys[i], val);
                if (vals[i] == null) {
                    vals[i] = new HashSet<>();
                }
                vals[i].add(val);
                assertEquals(1, hm.getCardinality(keys[i], val, val+1));
            }
        }

        for (int i = 0; i < n_tests; i++) {
            for (int j = 0; j < 1000; j++) {
                int x = rand.nextInt(130);
                int y = rand.nextInt(130);
                int from = Math.min(x, y);
                int to = Math.max(x, y);

                int cnt = 0;
                for (long k: vals[i]) {
                    if (from <= k && k < to) {
                        cnt++;
                    }
                }

                assertEquals(cnt, hm.getCardinality(keys[i], from, to));
            }
        }

        hm.reset();

        for (int i = 0; i < n_tests; i++) {
            assertEquals(0, hm.getCardinality(keys[i], 1, 100));
        }
    }


    private static int getPositionInt(Long2BitLongaHashMap.MapData curData, long key) {
        long h = HashCommon.murmurHash3(key);
        int pos = (int)(h & (long)curData.capacityMask);
        long[] keys = curData.keys;

        while(keys[pos] != 0L && keys[pos] != key) {
            ++pos;
            if (pos == keys.length) {
                pos = 0;
            }
        }

        return pos;
    }
}