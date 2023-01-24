package structures.map;

public class MutableLongBitShortaEntry {
        private long key;
        private short[] value;


        public long getKey() { return key; }
        public short[] getValue() { return value; }

        public void setKey(long key) { this.key = key; }
        public void setValue(short[] value) { this.value = value; }
}