package algo;

import java.util.ArrayList;
import java.util.List;

public class SingleNode {
    public String sequence;
    public int id;
    public boolean deleted;
    public SingleNode rc;
    public List<SingleNode> neighbors;

    public SingleNode(String sequence, int id) {
        this.sequence = sequence;
        this.id = id;

        this.deleted = false;
        this.neighbors = new ArrayList<>();
    }

}
