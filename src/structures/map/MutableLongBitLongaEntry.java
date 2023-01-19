package structures.map;

public class MutableLongBitLongaEntry {
        private long key;
        private long[] value;


        public long getKey() { return key; }
        public long[] getValue() { return value; }

        public void setKey(long key) { this.key = key; }
        public void setValue(long[] value) { this.value = value; }
}