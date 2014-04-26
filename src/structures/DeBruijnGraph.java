package structures;

import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;

/**
 * Vladimir Ulyantsev
 * Date: 05.04.14
 * Time: 0:09
 */
public class DeBruijnGraph {

    private int k;

    public ArrayLong2IntHashMap hm;

    private DeBruijnGraph(int k, ArrayLong2IntHashMap hm) {
        this.k = k;
        this.hm = hm;
    }

    public int getK() {
        return k;
    }
}
