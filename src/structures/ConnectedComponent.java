package structures;

import it.unimi.dsi.fastutil.longs.LongArrayList;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class ConnectedComponent implements Comparable<ConnectedComponent> {

    public long size;
    public long weight;


    // frequently used (not always)
    public int no;
    public int usedFreqThreshold;

    /**
     * Stores kmers if the component is small (less than b2 vertices).
     */
    public List<Long> kmers;

    /**
     * Is used in ComponentsBuilder if the component is big.
     */
    public BigLong2ShortHashMap hm = null;



    public ConnectedComponent() {
        kmers = new LongArrayList();
        size = 0;
        weight = 0;
    }


    public void add(long kmer, short w) {
        kmers.add(kmer);
        size++;
        weight += w;
    }
    public void add(long kmer) {
        kmers.add(kmer);
        size++;
    }



    public static void saveComponents(Collection<ConnectedComponent> components, String fp) throws IOException {
        DataOutputStream outputStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fp)));
        outputStream.writeInt(components.size());

        for (ConnectedComponent component : components) {
            outputStream.writeInt((int) component.size);
            outputStream.writeLong(component.weight);
            for (long kmer : component.kmers) {
                outputStream.writeLong(kmer);
            }
        }

        outputStream.close();
    }

    public static List<ConnectedComponent> loadComponents(File file) throws ExecutionFailedException {
        try {
            DataInputStream inputStream = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));
            int cnt = inputStream.readInt();
            List<ConnectedComponent> res = new ArrayList<ConnectedComponent>(cnt);

            for (int i = 0; i < cnt; i++) {
                int componentSize = inputStream.readInt();
                ConnectedComponent component = new ConnectedComponent();

                component.weight = inputStream.readLong();
                for (int j = 0; j < componentSize; j++) {
                    component.add(inputStream.readLong());
                }
                res.add(component);
                component.no = i + 1;
            }
            inputStream.close();
            return res;
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException("Can't load components: file not found", e);
        } catch (EOFException e) {
            throw new ExecutionFailedException("Can't load components: file corrupted or format mismatch! " +
                    "Do you set a wrong file?", e);
        } catch (IOException e) {
            throw new ExecutionFailedException("Can't load components: unknown IOException", e);
        }
    }


    @Override
    public int compareTo(ConnectedComponent o) {
        int sign = usedFreqThreshold - o.usedFreqThreshold;
        if (sign != 0) {
            return sign;
        }
        sign = -Long.compare(weight, o.weight);
        if (sign != 0) {
            return sign;
        }
        return -Long.compare(size, o.size);
    }
}
