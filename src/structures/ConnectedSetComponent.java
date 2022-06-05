package structures;

import org.apache.commons.lang.mutable.MutableLong;
import ru.ifmo.genetics.structures.set.LongHashSet;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Collection;

public class ConnectedSetComponent extends ConnectedComponent {
    public LongHashSet kmersSet;

    public ConnectedSetComponent() {
        kmersSet = new LongHashSet();
        size = 0;
        weight = 0;
    }

    public static void saveSetComponents(Collection<ConnectedSetComponent> components, String fp) throws IOException {
        DataOutputStream outputStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fp)));
        outputStream.writeInt(components.size());

        for (ConnectedSetComponent component : components) {
            outputStream.writeInt((int) component.size);
            outputStream.writeLong(component.weight);
            for (MutableLong kmer : component.kmersSet) {
                outputStream.writeLong(kmer.toLong());
            }
        }

        outputStream.close();
    }

    public boolean contains(long kmer) {
        return kmersSet.contains(kmer);
    }

    @Override
    public void add(long kmer, short w) {
        kmersSet.add(kmer);
        size++;
        weight += w;
    }

    @Override
    public void add(long kmer) {
        kmersSet.add(kmer);
        size++;
    }
}
