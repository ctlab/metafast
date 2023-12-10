package algo;

import java.io.Serializable;
import java.util.*;

class Observation {
    private final double value;
    private final int group;

    public Observation(double value, int group) {
        this.value = value;
        this.group = group;
    }

    public double getValue() {
        return this.value;
    }

    public int getGroup() {
        return this.group;
    }
}

class RankedObservation extends Observation {
    private double rank;

    public RankedObservation(double value, int group) {
        super(value, group);
    }

    public double getRank() {
        return this.rank;
    }

    public void setRank(double rank) {
        this.rank = rank;
    }
}

class ObservationComparator implements Comparator<RankedObservation>, Serializable {

    ObservationComparator() {
    }

    public int compare(RankedObservation o1, RankedObservation o2) {
        if (o1.getValue() < o2.getValue()) {
            return -1;
        } else {
            return o1.getValue() > o2.getValue() ? 1 : 0;
        }
    }
}

public class KruskalWallisTest {
    int numberOfGroups;
    Comparator<RankedObservation> comparator;
    List<RankedObservation> data;

    public KruskalWallisTest(int numberOfGroups) {
        this.numberOfGroups = numberOfGroups;
        this.comparator = new ObservationComparator();
        this.data = new ArrayList();
        if (numberOfGroups <= 1) {
            throw new IllegalArgumentException("requires two or more groups");
        }
    }

    public void add(double value, int group) {
        if (group >= 0 && group < this.numberOfGroups) {
            this.data.add(new RankedObservation(value, group));
        } else {
            throw new IllegalArgumentException();
        }

    }

    public void addAll(double[] values, int group) {
        double[] arr$ = values;
        int len$ = values.length;

        for (int i$ = 0; i$ < len$; ++i$) {
            double value = arr$[i$];
            this.add(value, group);
        }
    }

    double H() {
        int[] n = new int[this.numberOfGroups];
        double[] rbar = new double[this.numberOfGroups];

        int var10001;
        RankedObservation observation;
        for (Iterator i$ = this.data.iterator(); i$.hasNext(); rbar[var10001] += observation.getRank()) {
            observation = (RankedObservation) i$.next();
            ++n[observation.getGroup()];
            var10001 = observation.getGroup();
        }

        double H = 0.0D;

        int N;
        for (N = 0; N < this.numberOfGroups; ++N) {
            H += Math.pow(rbar[N], 2.0D) / (double) n[N];
        }

        N = this.data.size();
        return 12.0D / (double) (N * (N + 1)) * H - 3.0D * (double) (N + 1);
    }

    double C() {
        int N = this.data.size();
        double C = 0.0D;

        int j;
        for (int i = 0; i < N; i = j) {
            for (j = i + 1; j < N && ((RankedObservation) this.data.get(i)).getValue() == ((RankedObservation) this.data.get(j)).getValue(); ++j) {
            }

            C += Math.pow((double) (j - i), 3.0D) - (double) (j - i);
        }

        return 1.0D - C / (Math.pow((double) N, 3.0D) - (double) N);
    }

    public void update() {
        Collections.sort(this.data, this.comparator);

        int j;
        for (int i = 0; i < this.data.size(); i = j) {
            j = i + 1;

            double rank;
            for (rank = (double) (i + 1); j < this.data.size() && ((RankedObservation) this.data.get(i)).getValue() == ((RankedObservation) this.data.get(j)).getValue(); ++j) {
                rank += (double) (j + 1);
            }

            rank /= (double) (j - i);

            for (int k = i; k < j; ++k) {
                ((RankedObservation) this.data.get(k)).setRank(rank);
            }
        }

    }

    public double test() {
        this.update();
        double H = this.H();
        double C = this.C();
        if (C == 0.0D) {
            return 1;
        } else {
            return H / C;
        }
    }
}