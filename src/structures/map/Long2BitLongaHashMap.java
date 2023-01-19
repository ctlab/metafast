package structures.map;

import ru.ifmo.genetics.structures.set.LongHashSet;
import ru.ifmo.genetics.utils.NumUtils;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class Long2BitLongaHashMap extends LongHashSet implements Long2BitLongaHashMapInterface {

    private final static int BITS_PER_WORD = 6; // size(long)=64=2^6



    protected class MapData extends SetData {
        protected final long[][] values;
        protected volatile long[] valueForFreeKey;
        protected volatile int sizeBitSet;

        public MapData(int capacity, float maxLoadFactor, int sizeBitSet) {
            super(capacity, maxLoadFactor);
            values = new long[capacity][];
            /*for (int i=0; i<capacity; i++) {
                values[i] = new long[(sizeBitSet>>BITS_PER_WORD) + 1];
            }*/
            valueForFreeKey = new long[(sizeBitSet>>BITS_PER_WORD) + 1];
            this.sizeBitSet = sizeBitSet;
        }
    }

    protected volatile MapData data;

    // constructors
    public Long2BitLongaHashMap() {
        this(20, DEFAULT_MAX_LOAD_FACTOR, 0);  // 1 M elements
    }
    public Long2BitLongaHashMap(int capacity) {
        this(
                NumUtils.getPowerOf2(capacity),
                DEFAULT_MAX_LOAD_FACTOR,
                0
        );
    }
    public Long2BitLongaHashMap(int logCapacity, float maxLoadFactor, int sizeBitSet) {
        if (logCapacity > 30) {
            throw new IllegalArgumentException("log capacity > 30!");
        }

        this.maxLoadFactor = maxLoadFactor;
        int capacity = 1 << logCapacity;
        data = new MapData(capacity, maxLoadFactor, sizeBitSet);
        super.data = data;
    }

    @Override
    public boolean add(long key) {
        return set(key, 0) == null;
    }


    @Override
    public long[] set(long key, int bitIndex) {
        if (key == FREE) {
            writeLock.lock();
            try {
                MapData curData = data;
                long[] prev = curData.valueForFreeKey.clone();

                curData.valueForFreeKey[bitIndex>>BITS_PER_WORD]|=1L<<(bitIndex&((1L<<BITS_PER_WORD) - 1));
                if (!curData.containsFreeKey) {
                    curData.containsFreeKey = true;
                    curData.size++;
                    return null;
                }
                return prev;
            } finally {
                writeLock.unlock();
            }
        }

        while (true) {
            MapData curData = data;
            int pos = getPositionInt(curData, key);
            writeLock.lock();
            try {
                if (curData == data && (curData.keys[pos] == FREE || curData.keys[pos] == key)) {  // i.e. nothing has changed
                    if (curData.values[pos] == null) {
                        curData.values[pos] = new long[(data.sizeBitSet>>BITS_PER_WORD) + 1];
                    }
                    long[] prev = curData.values[pos].clone();
                    curData.values[pos][bitIndex>>BITS_PER_WORD]|=1L<<(bitIndex&((1L<<BITS_PER_WORD) - 1));
                    if (curData.keys[pos] == FREE) {
                        curData.keys[pos] = key;
                        curData.size++;
                        if (curData.size >= curData.maxFill) {
                            enlargeAndRehash();
                        }
                        return null;
                    }
                    return prev;
                }
            } finally {
                writeLock.unlock();
            }
        }
    }

    private void enlargeAndRehash() {
        MapData curData = data;
        if (curData.capacity > Integer.MAX_VALUE / 2) {
            throw new RuntimeException("Can't enlarge map (can't create single array of 2^31 elements)!");
        }
        int newCapacity = 2 * curData.capacity;
        MapData newData = new MapData(newCapacity, maxLoadFactor, curData.sizeBitSet);

        // coping elements
        for (int oldPos = 0; oldPos < curData.keys.length; oldPos++) {
            long key = curData.keys[oldPos];
            if (key != FREE) {
                int pos = getPositionInt(newData, key);
                newData.keys[pos] = key;
                newData.values[pos] = curData.values[oldPos];
            }
        }
        newData.containsFreeKey = curData.containsFreeKey;
        newData.valueForFreeKey = curData.valueForFreeKey;
        newData.size = curData.size;

        data = newData;
        super.data = newData;
    }

    @Override
    public boolean get(long key, int bitIndex) {
        long[] val = get(key);
        if (val == null) {
            return false;
        }
        return ((val[bitIndex>>BITS_PER_WORD]>>(bitIndex&((1L<<BITS_PER_WORD) - 1)))&1) == 1;
    }

    @Override
    public long[] get(long key) {
        MapData curData = data;
        if (key == FREE) {
            if (!curData.containsFreeKey) {
                return null;
            }
            return curData.valueForFreeKey;
        }
        int pos = getPositionInt(curData, key);
        if (curData.keys[pos] == key) {
            return curData.values[pos];
        } else {
            // assuming keys[pos] == FREE
            return null;
        }
    }

    @Override
    public long[] getWithEmpty(long key) {
        long[] value = get(key);
        if (value == null)
            return new long[1];
        return value;
    }

    @Override
    public int getCardinality(long key) {
        long[] value = get(key);
        if (value == null)
            return 0;
        int count = 0;
        for (long l : value) {
            count += countSetBits(l);
        }
        return count;
    }

    @Override
    public int getCardinality(long key, int from, int to) {
        long[] value = get(key);
        if (value == null)
            return 0;
        int count = 0;

        long bitFrom = from&((1L<<BITS_PER_WORD) - 1);
        long bitTo = to&((1L<<BITS_PER_WORD) - 1);
        if ((from>>BITS_PER_WORD) == (to>>BITS_PER_WORD)) {
            return countSetBits(value[from>>BITS_PER_WORD] & ((1L<<bitTo)-(1L<<bitFrom)));
        } else {
            long all = 0xFFFFFFFFFFFFFFFFL;
            count += countSetBits(value[from>>BITS_PER_WORD] & (all - ((1L<<bitFrom)-1)));
            for (int i = (from>>BITS_PER_WORD) + 1; i < (to>>BITS_PER_WORD); i++) {
                count += countSetBits(value[i]);
            }
            count += countSetBits(value[to>>BITS_PER_WORD] & ((1L<<bitTo) - 1));
        }
        return count;
    }

    private int countSetBits(long i)
    {
        i = i - ((i >> 1) & 0x5555555555555555L);
        i = (i & 0x3333333333333333L) +
                ((i >> 2) & 0x3333333333333333L);
        i = ((i + (i >> 4)) & 0x0F0F0F0F0F0F0F0FL);
        return (int)((i*0x0101010101010101L)>>56);
    }

    /*
    Inefficient implementation
    @Override
    public int getCardinality(long key) {
        long[] value = get(key);
        if (value == null)
            return 0;
        return getCardinality(value, 0, value.length<<BITS_PER_WORD);
    }



    @Override
    public int getCardinality(long key, int from, int to) {
        long[] value = get(key);
        if (value == null)
            return 0;
        return getCardinality(value, from, to);
    }

    private int getCardinality(long[] bits, int from, int to) {
        int counter = 0;
        for (int i = from; i < to; i++) {
            counter += (bits[i>>BITS_PER_WORD]>>(i&((1L<<BITS_PER_WORD) - 1)))&1;
        }
        return counter;
    }
    */

    @Override
    public boolean contains(long key) {
        return contains(data, key);
    }

    @Override
    public long size() { return data.size; }

    @Override
    public long capacity() { return data.capacity; }


    @Override
    public void reset() {
        writeLock.lock();
        MapData curData = data;
        try {
            Arrays.fill(curData.keys, FREE);
            Arrays.fill(curData.values, null);
            curData.containsFreeKey = false;
            curData.valueForFreeKey = new long[(curData.sizeBitSet>>BITS_PER_WORD) + 1];
            curData.size = 0;
        } finally {
            writeLock.unlock();
        }
    }

    @Override
    public void resetValues() {
        writeLock.lock();
        MapData curData = data;
        try {
            Arrays.fill(curData.values, null);
            curData.valueForFreeKey = new long[(curData.sizeBitSet>>BITS_PER_WORD) + 1];
        } finally {
            writeLock.unlock();
        }
    }


    @Override
    public long keyAt(long pos) {
        return elementAt(pos);
    }

    @Override
    public long[] valueAt(long pos) {
        MapData curData = data;
        if (pos == curData.capacity) {
            if (!curData.containsFreeKey) {
                return null;
            }
            return curData.valueForFreeKey;
        }

        if (curData.keys[(int) pos] == FREE) {
            return null;
        }
        return curData.values[(int) pos];
    }


    @Override
    public void write(DataOutput out) throws IOException {
        MapData curData = data;

        out.writeInt(curData.sizeBitSet);
        out.writeInt(curData.capacity);
        out.writeInt(curData.size);
        out.writeFloat(maxLoadFactor);

        for (int i = 0; i < curData.capacity; i++) {
            out.writeLong(curData.keys[i]);
            out.writeInt(curData.values[i].length);
            for (int j = 0; j < curData.values[i].length; j++) {
                out.writeLong(curData.values[i][j]);
            }
        }
        out.writeBoolean(curData.containsFreeKey);
        out.writeInt(curData.valueForFreeKey.length);
        for (int j = 0; j < curData.valueForFreeKey.length; j++) {
            out.writeLong(curData.valueForFreeKey[j]);
        }
    }

    @Override
    public void readFields(DataInput in) throws IOException {
        int sizeBitSet = in.readInt();
        int capacity = in.readInt();
        int size = in.readInt();
        maxLoadFactor = in.readFloat();


        MapData newData = new MapData(capacity, maxLoadFactor, sizeBitSet);

        for (int i = 0; i < capacity; i++) {
            newData.keys[i] = in.readLong();
            int len = in.readInt();
            long[] bytes = new long[len];
            for (int j = 0; j < len; j++) {
                bytes[j] = in.readLong();
            }
            newData.values[i] = bytes;
        }
        newData.containsFreeKey = in.readBoolean();
        int len = in.readInt();
        long[] bytes = new long[len];
        for (int j = 0; j < len; j++) {
            bytes[j] = in.readLong();
        }
        newData.valueForFreeKey = bytes;
        newData.size = size;

        data = newData;
        super.data = newData;
    }

    @Override
    public Iterator<MutableLongBitLongaEntry> entryIterator() {
        return new MyIterator(data);
    }

    protected class MyIterator implements Iterator<MutableLongBitLongaEntry> {
        private final MapData curData;
        private int index = 0;
        private final MutableLongBitLongaEntry entry = new MutableLongBitLongaEntry();

        MyIterator(MapData curData) {
            this.curData = curData;
        }

        @Override
        public boolean hasNext() {
            while ((index < curData.capacity) && (curData.keys[index] == FREE)) {
                index++;
            }
            if (index < curData.capacity) {
                return true;
            }
            if (index == curData.capacity && curData.containsFreeKey) {
                return true;
            }
            return false;
        }

        @Override
        public MutableLongBitLongaEntry next() {
            if (hasNext()){
                if (index < curData.capacity) {
                    entry.setKey(curData.keys[index]);
                    entry.setValue(curData.values[index]);
                }
                if (index == curData.capacity) {
                    entry.setKey(FREE);
                    entry.setValue(curData.valueForFreeKey);
                }
                index++;
                return entry;
            }
            throw new NoSuchElementException();
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }
}