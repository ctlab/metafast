package algo;

import ru.ifmo.genetics.executors.BlockingThreadPoolExecutor;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import structures.Sequence;

import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

public class SequencesFinders {

    public static List<Sequence> thresholdStrategy(ArrayLong2IntHashMap hm,
                                                   int availableProcessors,
                                                   int freqThreshold,
                                                   int lenThreshold,
                                                   int k) throws InterruptedException {
        List<Sequence> ans = new CopyOnWriteArrayList<Sequence>();

        BlockingThreadPoolExecutor executor = new BlockingThreadPoolExecutor(availableProcessors);

        for (int i = 0; i < hm.hm.length; ++i) {
            executor.blockingExecute(new
                    AddSequencesShiftingRightTask(hm, hm.hm[i], freqThreshold, ans, lenThreshold, k));
        }

        System.out.println(executor.getTaskCount());
        executor.shutdownAndAwaitTermination();
        return ans;
    }
}
