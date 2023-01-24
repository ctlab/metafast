package structures.map;

import ru.ifmo.genetics.structures.set.LongHashSet;
import ru.ifmo.genetics.utils.NumUtils;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class Long2BitSetHashMap extends LongHashSet implements Long2BitSetHashMapInterface {

    protected class MapData extends SetData {
        protected final BitSet[] values;
        protected volatile BitSet valueForFreeKey;
        protected volatile int sizeBitSet;

        public MapData(int capacity, float maxLoadFactor, int sizeBitSet) {
            super(capacity, maxLoadFactor);
            values = new BitSet[capacity];
            for (int i=0; i<capacity; i++) {
                values[i] = new BitSet();
            }
            //valueForFreeKey = new BitSet(sizeBitSet);
            valueForFreeKey = new BitSet();
            this.sizeBitSet = sizeBitSet;
        }
    }

    protected volatile Long2BitSetHashMap.MapData data;

    // constructors
    public Long2BitSetHashMap() {
        this(20, DEFAULT_MAX_LOAD_FACTOR, 0);  // 1 M elements
    }
    public Long2BitSetHashMap(int capacity) {
        this(
                NumUtils.getPowerOf2(capacity),
                DEFAULT_MAX_LOAD_FACTOR,
                0
        );
    }
    public Long2BitSetHashMap(int logCapacity, float maxLoadFactor, int sizeBitSet) {
        if (logCapacity > 30) {
            throw new IllegalArgumentException("log capacity > 30!");
        }

        this.maxLoadFactor = maxLoadFactor;
        int capacity = 1 << logCapacity;
        data = new Long2BitSetHashMap.MapData(capacity, maxLoadFactor, sizeBitSet);
        super.data = data;
    }

    @Override
    public boolean add(long key) {
        return set(key, 0) == null;
    }


    @Override
    public BitSet set(long key, int bitIndex) {
        if (key == FREE) {
            writeLock.lock();
            try {
                MapData curData = data;
                BitSet prev = curData.valueForFreeKey;
                curData.valueForFreeKey.set(bitIndex);
                if (!curData.containsFreeKey) {
                    curData.containsFreeKey = true;
                    curData.size++;
                }
                return prev;
            } finally {
                writeLock.unlock();
            }
        }

        while (true) {
            Long2BitSetHashMap.MapData curData = data;
            int pos = getPositionInt(curData, key);
            writeLock.lock();
            try {
                if (curData == data && (curData.keys[pos] == FREE || curData.keys[pos] == key)) {  // i.e. nothing has changed
                    BitSet prev = curData.values[pos];
                    curData.values[pos].set(bitIndex);
                    if (curData.keys[pos] == FREE) {
                        curData.keys[pos] = key;
                        curData.size++;
                        if (curData.size >= curData.maxFill) {
                            enlargeAndRehash();
                        }
                    }
                    return prev;
                }
            } finally {
                writeLock.unlock();
            }
        }
    }

    private void enlargeAndRehash() {
        Long2BitSetHashMap.MapData curData = data;
        if (curData.capacity > Integer.MAX_VALUE / 2) {
            throw new RuntimeException("Can't enlarge map (can't create single array of 2^31 elements)!");
        }
        int newCapacity = 2 * curData.capacity;
        Long2BitSetHashMap.MapData newData = new Long2BitSetHashMap.MapData(newCapacity, maxLoadFactor, curData.sizeBitSet);

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
        return getWithEmpty(key).get(bitIndex);
    }

    @Override
    public BitSet get(long key) {
        Long2BitSetHashMap.MapData curData = data;
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
    public BitSet getWithEmpty(long key) {
        BitSet value = get(key);
        if (value == null)
            return new BitSet();
        return value;
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
        Long2BitSetHashMap.MapData curData = data;
        try {
            Arrays.fill(curData.keys, FREE);
            Arrays.fill(curData.values, null);
            curData.containsFreeKey = false;
            //curData.valueForFreeKey = new BitSet(curData.sizeBitSet);
            curData.valueForFreeKey = new BitSet();
            curData.size = 0;
        } finally {
            writeLock.unlock();
        }
    }

    @Override
    public void resetValues() {
        writeLock.lock();
        Long2BitSetHashMap.MapData curData = data;
        try {
            Arrays.fill(curData.values, null);
            //curData.valueForFreeKey = new BitSet(curData.sizeBitSet);
            curData.valueForFreeKey = new BitSet();
        } finally {
            writeLock.unlock();
        }
    }


    @Override
    public long keyAt(long pos) {
        return elementAt(pos);
    }

    @Override
    public BitSet valueAt(long pos) {
        Long2BitSetHashMap.MapData curData = data;
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
        Long2BitSetHashMap.MapData curData = data;

        out.writeInt(curData.sizeBitSet);
        out.writeInt(curData.capacity);
        out.writeInt(curData.size);
        out.writeFloat(maxLoadFactor);

        for (int i = 0; i < curData.capacity; i++) {
            out.writeLong(curData.keys[i]);
            out.writeInt(curData.values[i].length());
            out.write(curData.values[i].toByteArray());
        }
        out.writeBoolean(curData.containsFreeKey);
        out.writeInt(curData.valueForFreeKey.length());
        out.write(curData.valueForFreeKey.toByteArray());
    }

    @Override
    public void readFields(DataInput in) throws IOException {
        int sizeBitSet = in.readInt();
        int capacity = in.readInt();
        int size = in.readInt();
        maxLoadFactor = in.readFloat();


        Long2BitSetHashMap.MapData newData = new Long2BitSetHashMap.MapData(capacity, maxLoadFactor, sizeBitSet);

        for (int i = 0; i < capacity; i++) {
            newData.keys[i] = in.readLong();
            int len = in.readInt();
            byte[] bytes = new byte[len];
            in.readFully(bytes);
            newData.values[i] = BitSet.valueOf(bytes);
        }
        newData.containsFreeKey = in.readBoolean();
        int len = in.readInt();
        byte[] bytes = new byte[len];
        in.readFully(bytes);
        newData.valueForFreeKey = BitSet.valueOf(bytes);
        newData.size = size;

        data = newData;
        super.data = newData;
    }

    @Override
    public Iterator<MutableLongBitSetEntry> entryIterator() {
        return new Long2BitSetHashMap.MyIterator(data);
    }

    protected class MyIterator implements Iterator<MutableLongBitSetEntry> {
        private final Long2BitSetHashMap.MapData curData;
        private int index = 0;
        private final MutableLongBitSetEntry entry = new MutableLongBitSetEntry();

        MyIterator(Long2BitSetHashMap.MapData curData) {
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
        public MutableLongBitSetEntry next() {
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