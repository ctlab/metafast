package io;

import ru.ifmo.genetics.dna.Dna;

import java.util.List;
import java.util.concurrent.CountDownLatch;

public abstract class ReadsWorker implements Runnable {

    private ReadsDispatcher dispatcher = null;
    private CountDownLatch latch = null;

    boolean interrupted = false;


    void setDispatcher(ReadsDispatcher dispatcher) {
        this.dispatcher = dispatcher;
    }
    void setLatch(CountDownLatch latch) {
        this.latch = latch;
    }



    public abstract void process(List<Dna> reads);


    @Override
    public void run() {
        if (dispatcher == null || latch == null) {
            throw new RuntimeException("Not full initialization!");
        }
        while (!interrupted) {
            List<Dna> list = dispatcher.getWorkRange();
            if (list == null) {
                break;
            }
            process(list);
        }
        latch.countDown();
    }

    public void interrupt() {
        interrupted = true;
    }
}
