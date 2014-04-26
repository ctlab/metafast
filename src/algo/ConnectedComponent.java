package algo;

import it.unimi.dsi.fastutil.longs.LongArrayList;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class ConnectedComponent {
    public static void saveComponents(Collection<ConnectedComponent> components, String fp) throws IOException {
        DataOutputStream outputStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fp)));
        outputStream.writeInt(components.size());

        for (ConnectedComponent component : components) {
            outputStream.writeInt(component.size());
            outputStream.writeLong(component.getWeight());
            for (long kmer : component.kmers) {
                outputStream.writeLong(kmer);
            }
        }

        outputStream.close();
    }

    public static List<ConnectedComponent>
    loadComponents(String fp) throws IOException {
        DataInputStream inputStream = new DataInputStream(new BufferedInputStream(new FileInputStream(fp)));
        int componentsCnt = inputStream.readInt();
        List<ConnectedComponent> res = new ArrayList<ConnectedComponent>(componentsCnt);

        for (int i = 0; i < componentsCnt; i++) {
            int componentSize = inputStream.readInt();
            ConnectedComponent component = new ConnectedComponent();

            component.setWeight(inputStream.readLong());
            for (int j = 0; j < componentSize; j++) {
                component.add(inputStream.readLong());
            }
            res.add(component);
        }
        inputStream.close();

        return res;
    }

    public List<Long> kmers;

    private long weight;

    public ConnectedComponent() {
        kmers = new LongArrayList();
    }

    public int size() {
        return kmers.size();
    }

    public void add(long kmer) {
        kmers.add(kmer);
    }

    public void setWeight(long weight) {
        this.weight = weight;
    }

    public long getWeight() {
        return weight;
    }
}
