package structures;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;

//todo  убрать копипасту
public class ConnectedSetComponent extends ConnectedComponent{
    public HashSet<Long> kmersset;
    public boolean contains(long kmer) {
        return kmersset.contains(kmer);
    }

    public ConnectedSetComponent() {
        kmersset = new HashSet<>();
        size = 0;
        weight = 0;
    }

    public ConnectedSetComponent(SequenceComponent component) {
        kmersset = new HashSet<>();
        kmersset.addAll(component.kmers);
        size = component.size;
        weight = component.weight;
    }

    @Override
    public void add(long kmer, short w) {
        kmersset.add(kmer);
        size++;
        weight += w;
    }

    @Override
    public void add(long kmer) {
        kmersset.add(kmer);
        size++;
    }


    public static void saveSetComponents(Collection<ConnectedSetComponent> components, String fp) throws IOException {
        DataOutputStream outputStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fp)));
        outputStream.writeInt(components.size());

        for (ConnectedSetComponent component : components) {
            outputStream.writeInt((int) component.size);
            outputStream.writeLong(component.weight);
            for (long kmer : component.kmersset) {
                outputStream.writeLong(kmer);
            }
        }

        outputStream.close();
    }
}
