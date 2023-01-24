package structures.map;

import org.apache.commons.lang.mutable.MutableLong;
import org.apache.hadoop.io.Writable;

import java.util.BitSet;
import java.util.Iterator;

public interface Long2BitSetHashMapInterface extends Writable, Iterable<MutableLong> {
    /**
     * @return the previous value of this key, or null, if no value was associated.
     */
    public BitSet set(long key, int bitIndex);

    /**
     * @return the value of the bit with the specified @bitIndex
     */
    public boolean get(long key, int bitIndex);

    /**
     * @return null, if not found
     */
    public BitSet get(long key);

    /**
     * @return empty BitSet, if not found
     */
    public BitSet getWithEmpty(long key);
    public boolean contains(long key);


    public long size();
    public long capacity();

    public void reset();
    public void resetValues();

    public Iterator<MutableLongBitSetEntry> entryIterator();


    // Methods, that use information about the internal structure. May be unsupported.
    /**
     * Call this method to prepare to the future requests.
     * Assuming no other thread modifying map!
     */
    public void prepare();
    public long maxPosition();
    /**
     * @return pos or -1, if not found
     */
    public long getPosition(long key);
    /**
     * Returns FREE, if pos is free.
     */
    public long keyAt(long pos);
    /**
     * Returns null, if pos is free.
     */
    public BitSet valueAt(long pos);
    public boolean containsAt(long pos);

}