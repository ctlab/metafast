package structures;

import org.apache.commons.lang.mutable.MutableLong;
import ru.ifmo.genetics.structures.set.LongHashSet;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;

//todo  убрать копипасту
public class ConnectedSetComponent extends ConnectedComponent{
    public LongHashSet kmersset;
    public boolean contains(long kmer) {
        return kmersset.contains(kmer);
    }

    public ConnectedSetComponent() {
        kmersset = new LongHashSet();
        size = 0;
        weight = 0;
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
            for (MutableLong kmer : component.kmersset) {
                outputStream.writeLong(kmer.toLong());
            }
        }

        outputStream.close();
    }
}
