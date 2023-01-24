package structures.map;

import java.util.BitSet;

public class MutableLongBitSetEntry {

    private long key;
    private BitSet value;


    public long getKey() { return key; }
    public BitSet getValue() { return value; }

    public void setKey(long key) { this.key = key; }
    public void setValue(BitSet value) { this.value = value; }
}
