package structures.map;

import ru.ifmo.genetics.structures.set.LongHashSet;
import ru.ifmo.genetics.utils.NumUtils;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class Long2BitShortaHashMap extends LongHashSet implements Long2BitShortaHashMapInterface {

    private final static int BITS_PER_WORD = 4; // size(short)=16=2^4



    protected class MapData extends SetData {
        protected final short[][] values;
        protected volatile short[] valueForFreeKey;
        protected volatile int sizeBitSet;

        public MapData(int capacity, float maxLoadFactor, int sizeBitSet) {
            super(capacity, maxLoadFactor);
            values = new short[capacity][];
            valueForFreeKey = new short[(sizeBitSet>>BITS_PER_WORD) + 1];
            this.sizeBitSet = sizeBitSet;
        }
    }

    protected volatile MapData data;

    // constructors
    public Long2BitShortaHashMap() {
        this(20, DEFAULT_MAX_LOAD_FACTOR, 0);  // 1 M elements
    }
    public Long2BitShortaHashMap(int capacity) {
        this(
                NumUtils.getPowerOf2(capacity),
                DEFAULT_MAX_LOAD_FACTOR,
                0
        );
    }
    public Long2BitShortaHashMap(int logCapacity, float maxLoadFactor, int sizeBitSet) {
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
    public short[] set(long key, int bitIndex) {
        if (key == FREE) {
            writeLock.lock();
            try {
                MapData curData = data;
                short[] prev = curData.valueForFreeKey.clone();

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
                        curData.values[pos] = new short[(data.sizeBitSet>>BITS_PER_WORD) + 1];
                    }
                    short[] prev = curData.values[pos].clone();
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

        // copying elements
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
        short[] val = get(key);
        if (val == null) {
            return false;
        }
        return ((val[bitIndex>>BITS_PER_WORD]>>(bitIndex&((1L<<BITS_PER_WORD) - 1)))&1) == 1;
    }

    @Override
    public short[] get(long key) {
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
    public short[] getWithEmpty(long key) {
        short[] value = get(key);
        if (value == null)
            return new short[1];
        return value;
    }

    @Override
    public int getCardinality(long key) {
        short[] value = get(key);
        if (value == null)
            return 0;
        int count = 0;
        for (short l : value) {
            count += countSetBits(l);
        }
        return count;
    }

    @Override
    public int getCardinality(long key, int from, int to) {
        short[] value = get(key);
        if (value == null)
            return 0;
        int count = 0;

        long bitFrom = from&((1L<<BITS_PER_WORD) - 1);
        long bitTo = to&((1L<<BITS_PER_WORD) - 1);
        if ((from>>BITS_PER_WORD) == (to>>BITS_PER_WORD)) {
            return countSetBits((short) (value[from>>BITS_PER_WORD] & ((1L<<bitTo)-(1L<<bitFrom))));
        } else {
            short all = (short)0xFFFF;
            count += countSetBits((short) (value[from>>BITS_PER_WORD] & (all - ((1L<<bitFrom)-1))));
            for (int i = (from>>BITS_PER_WORD) + 1; i < (to>>BITS_PER_WORD); i++) {
                count += countSetBits(value[i]);
            }
            count += countSetBits((short) (value[to>>BITS_PER_WORD] & ((1L<<bitTo) - 1)));
        }
        return count;
    }

    private int countSetBits(short i) {
        i = (short) ((i & 0x5555) + ((i >> 1) & 0x5555));
        i = (short) ((i & 0x3333) + ((i >> 2) & 0x3333));
        i = (short) ((i & 0x0F0F) + ((i >> 4) & 0x0F0F));
        return ((short)(i * 0x0101))>>8;
    }

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
            curData.valueForFreeKey = new short[(curData.sizeBitSet>>BITS_PER_WORD) + 1];
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
            curData.valueForFreeKey = new short[(curData.sizeBitSet>>BITS_PER_WORD) + 1];
        } finally {
            writeLock.unlock();
        }
    }


    @Override
    public long keyAt(long pos) {
        return elementAt(pos);
    }

    @Override
    public short[] valueAt(long pos) {
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
                out.writeShort(curData.values[i][j]);
            }
        }
        out.writeBoolean(curData.containsFreeKey);
        out.writeInt(curData.valueForFreeKey.length);
        for (int j = 0; j < curData.valueForFreeKey.length; j++) {
            out.writeShort(curData.valueForFreeKey[j]);
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
            short[] bytes = new short[len];
            for (int j = 0; j < len; j++) {
                bytes[j] = in.readShort();
            }
            newData.values[i] = bytes;
        }
        newData.containsFreeKey = in.readBoolean();
        int len = in.readInt();
        short[] bytes = new short[len];
        for (int j = 0; j < len; j++) {
            bytes[j] = in.readShort();
        }
        newData.valueForFreeKey = bytes;
        newData.size = size;

        data = newData;
        super.data = newData;
    }

    @Override
    public Iterator<MutableLongBitShortaEntry> entryIterator() {
        return new MyIterator(data);
    }

    protected class MyIterator implements Iterator<MutableLongBitShortaEntry> {
        private final MapData curData;
        private int index = 0;
        private final MutableLongBitShortaEntry entry = new MutableLongBitShortaEntry();

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
        public MutableLongBitShortaEntry next() {
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