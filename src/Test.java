import algo.KmerOperations;
import it.unimi.dsi.fastutil.longs.LongArrayFIFOQueue;
import ru.ifmo.genetics.dna.kmers.ShortKmer;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class Test {

    public static void main(String[] args) {
        try {
            new FileReader(new File("xx"));
        } catch (IOException e) {
            System.err.println("Error while reading: \n" + e);
        }
    }


}
