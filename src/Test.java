import algo.KmerOperations;
import it.unimi.dsi.fastutil.longs.LongArrayFIFOQueue;
import ru.ifmo.genetics.dna.kmers.ShortKmer;

/**
 * Vladimir Ulyantsev
 * Date: 19.02.14
 * Time: 14:02
 */
public class Test {

    public static void main(String[] args) {
        ShortKmer kmer = new ShortKmer("ATAGC");
        for (long neighbour : KmerOperations.possibleNeighbours(kmer.fwKmer(), 5)) {
            System.out.print(new ShortKmer(neighbour, 5).toString() + " ");
        }
        System.out.println();

        kmer = new ShortKmer("GCTAT");
        for (long neighbour : KmerOperations.possibleNeighbours(kmer.fwKmer(), 5)) {
            System.out.print(new ShortKmer(neighbour, 5).toString() + " ");
        }
        System.out.println();

        LongArrayFIFOQueue queue = new LongArrayFIFOQueue();
        queue.enqueue(1L);
        queue.enqueue(239);
        System.out.println(queue.dequeue());
        System.out.println(queue.dequeue());
    }


}
