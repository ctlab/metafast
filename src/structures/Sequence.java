package structures;


import ru.ifmo.genetics.dna.kmers.ShortKmer;

public class Sequence {
    String repr;

    long totalWeight, minWeight, maxWeight;

    public Sequence(String repr, long totalWeight, long minWeight, long maxWeight) {
        this.repr = repr;
        this.totalWeight = totalWeight;
        this.minWeight = minWeight;
        this.maxWeight = maxWeight;
    }

    public int length() {
        return repr.length();
    }

    @Override
    public String toString() {
        return repr;
    }

    public int averageWeight() {
        return (int) (totalWeight / length());
    }

    public ShortKmer startKmer(int k) {
        return new ShortKmer(repr.substring(0, k));
    }

    public ShortKmer endKmer(int k) {
        return new ShortKmer(repr.substring(repr.length() - k, repr.length()));
    }
}
