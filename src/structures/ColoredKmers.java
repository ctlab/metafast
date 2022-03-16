package structures;

import it.unimi.dsi.fastutil.longs.LongArrayList;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;

import java.io.*;
import java.util.*;

public class ColoredKmers {


    private void add_color_int(long kmer, int color, int val) {
        long intaddval = (long) Math.pow(degree, color) + val;
        kmersColors.put(kmer, kmersColors.get(kmer) + intaddval);
    }

    public Integer[]  get_color_from_int(long kmer) {
        Integer[] res = new Integer[colorsCNT];
            long longres = kmersColors.get(kmer);
        for (int color = 0; color<colorsCNT; color++) {
            long curdegree = (long) Math.pow(degree, color) ;
//            long prevdegree = (long) Math.pow(degree, color + 1) ;
            int realv = (int) (longres / curdegree % degree);
            res[color] = realv;
        }
        return res;
    }

    public List<Long> kmers;
    public HashMap<Long, Long> kmersColors;
    //max cnt for color is 999999, max colorcnt now is 3
    public int colorsCNT;
    public long size;
    public long weight;
    private final int degree = 100000;
    private final double MIN_TO_COLOR = 0.75;

    public ColoredKmers(int colorsCNT) {
        this.colorsCNT = colorsCNT;
        kmers = new LongArrayList();
        size = 0;
        weight = 0;
        kmersColors = new HashMap<Long, Long>();
    }

    public void addColor(long kmer, int color) {
        addColor(kmer, color, 1);
    }

    public void addColor(long kmer, int color, int val) {
        if (!kmersColors.containsKey(kmer)) {
            kmers.add(kmer);
            size += 1;
            kmersColors.put(kmer, 0L);
        }
        add_color_int(kmer, color, val);
//        kmersColors.get(kmer)[color] += val;
    }


    private int argmax(Integer[] arr) {
        long sum = 0;
        int mi = 0;
        int mv = arr[0];
        for (int i = 0; i < arr.length; i++) {
            sum+=arr[i];
            if (arr[i] > mv) {
                mv = arr[i];
                mi = i;
            }
        }
        return (1.0 * mv / sum > MIN_TO_COLOR ) ? mi : -1;
    }

    public int getColor(long kmer) {
        int res = colorsCNT;
        if (kmersColors.containsKey(kmer)) {
            int vi = argmax(get_color_from_int(kmer));
            if (vi != -1) {
                res = vi;
            }
        }
        return res;
    }

    public Map<Long, Integer> getColors() {
        Map<Long, Integer> res = new HashMap<Long, Integer>();
        for (long kmer : kmers) {
            int v = getColor(kmer);
            res.put(kmer, v);
        }
        return res;
    }


    public void saveColorDouble(String fp) throws IOException {
        saveArrToFile(fp, false);
    }

    public void saveColorInt(String fp) throws IOException {
        saveArrToFile(fp, true);
    }

    public ColoredKmers(File file) throws ExecutionFailedException {
        try {
            DataInputStream inputStream = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));

            int size = inputStream.readInt();
            this.colorsCNT = inputStream.readInt();
            this.kmers = new ArrayList<Long>();
            this.kmersColors = new HashMap<Long, Long>();
            for (int j = 0; j < size; j++) {
                long kmer = inputStream.readLong();
                this.kmers.add(kmer);
                for (int k = 0; k < this.colorsCNT; k++) {
                    this.addColor(kmer, k, inputStream.readInt());
                }
            }
            System.out.println(this.size + " " + size);
            inputStream.close();
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException("Can't load components: file not found", e);
        } catch (EOFException e) {
            throw new ExecutionFailedException("Can't load components: file corrupted or format mismatch! " +
                    "Do you set a wrong file?", e);
        } catch (IOException e) {
            throw new ExecutionFailedException("Can't load components: unknown IOException", e);
        }
    }

    private void saveArrToFile(String fp, boolean norm) throws IOException {
        DataOutputStream outputStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fp)));
        outputStream.writeInt((int) this.size);
        outputStream.writeInt(this.colorsCNT);
        for (long kmer : this.kmers) {
            outputStream.writeLong(kmer);
            if (norm) {
                for (Integer v : get_color_from_int(kmer)) {
                    outputStream.writeInt(v);
                }
            } else {
                int normsum = 0;
                for (Integer v :  get_color_from_int(kmer)) {
                    normsum += v;
                }
                for (Integer v :  get_color_from_int(kmer)) {
                    outputStream.writeDouble(1.0 * v / normsum);
                }
            }
        }
        outputStream.close();
    }


}
