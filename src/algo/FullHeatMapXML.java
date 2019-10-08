package algo;

import com.sun.org.apache.xml.internal.serialize.Method;
import com.sun.org.apache.xml.internal.serialize.OutputFormat;
import com.sun.org.apache.xml.internal.serialize.XMLSerializer;
import org.apache.batik.anim.dom.SVGDOMImplementation;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import java.awt.*;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Constructing full heat map with dendrogram for n objects.
 */
public class FullHeatMapXML {

    // ===========================   Fixed params   ===========================

    final String svgNS = "http://www.w3.org/2000/svg";
    
    final int cellSize = 40;
    final int borderLineSize = 2;
    final int gridSize;

    final int dx_before_dendrogram = 100;
    final int dx_for_dendrogram;
    final int dx_for_names = 350;
    final int dy_for_names = 200;
    final int dy_before_color = 50;
    final int dy_after_color = 100;

    final int dx_scale_text = 120;
    final int dx_scale = 300;
    final int dy_scale = 30;


    public Color lowColor = Color.WHITE;
    public Color highColor = new Color(12, 61, 138);
    public boolean drawInnerLines = true;
    public Color innerLinesColor = Color.GRAY;

    final String fontXML = "bold 16px sans-serif";


    // ===========================   Variables   ==============================

    public final int n;
    public final String[] names;
    /**
     * Distance matrix with values from 0 to 1
     */
    protected final double[][] distMatrix;
    protected final double low, high;

    public final int[] perm;


    public FullHeatMapXML(double[][] distMatrix, String[] names) {
        this(distMatrix, 0, 1, names, false);
    }

    public FullHeatMapXML(double[][] distMatrix, String[] names, boolean invertColors) {
        this(distMatrix, 0, 1, names, invertColors);
    }

    public FullHeatMapXML(double[][] distMatrix, double low, double high, String[] names, boolean invertColors) {
        n = distMatrix.length;
        this.distMatrix = distMatrix;
        this.low = low;
        this.high = high;
        this.names = names;
        gridSize = n * cellSize + (n + 1) * borderLineSize;
        dx_for_dendrogram = getDendrogramSize(n);
        perm = new int[n];
        for (int i = 0; i < n; i++) {
            perm[i] = i;
        }
        if (invertColors) {
            Color tmp = highColor;
            highColor = lowColor;
            lowColor = tmp;
            innerLinesColor = Color.LIGHT_GRAY;
        }
    }


    static int getDendrogramSize(int n) {
        if (n < 2) {
            return 0;
        }
        int cn = 2;
        int cs = 50;
        while (cn < n) {
            cn++;
            if (3 <= cn && cn <= 5) {
                cs += 50;
            }
            if (6 <= cn && cn <= 10) {
                cs += 30;
            }
            if (11 <= cn && cn <= 20) {
                cs += 15;
            }
        }
        return cs;
    }


    // ==========================   Drawing heat map   ============================

    public Document createHeatMap() {
        return createHeatMap(0, 0);
    }

    public Document createHeatMap(int ddx, int ddy) {
        DOMImplementation impl = SVGDOMImplementation.getDOMImplementation();
        Document doc = impl.createDocument(svgNS, "svg", null);

        Element svgRoot = doc.getDocumentElement();
        svgRoot.setAttributeNS(null, "width", Integer.toString(gridSize));
        svgRoot.setAttributeNS(null, "height", Integer.toString(gridSize));

        { // background
            Element rectangle = doc.createElementNS(svgNS, "rect");
            rectangle.setAttributeNS(null, "x", Integer.toString(ddx));
            rectangle.setAttributeNS(null, "y", Integer.toString(ddy));
            rectangle.setAttributeNS(null, "width", Integer.toString(gridSize));
            rectangle.setAttributeNS(null, "height", Integer.toString(gridSize));
            rectangle.setAttributeNS(null, "fill", "white");
            svgRoot.appendChild(rectangle);
        }

        { // border lines
            int left = borderLineSize / 2;
            int right = n * (cellSize + borderLineSize) + borderLineSize / 2;

            Element border_1 = doc.createElementNS(svgNS, "line");
            border_1.setAttributeNS(null, "x1", Integer.toString(ddx));
            border_1.setAttributeNS(null, "y1", Integer.toString(left + ddy));
            border_1.setAttributeNS(null, "x2", Integer.toString(gridSize + ddx));
            border_1.setAttributeNS(null, "y2", Integer.toString(left + ddy));
            border_1.setAttributeNS(null, "stroke", "gray");
            border_1.setAttributeNS(null, "stroke-width", Integer.toString(borderLineSize));
            svgRoot.appendChild(border_1);

            Element border_2 = doc.createElementNS(svgNS, "line");
            border_2.setAttributeNS(null, "x1", Integer.toString(left + ddx));
            border_2.setAttributeNS(null, "y1", Integer.toString(ddy));
            border_2.setAttributeNS(null, "x2", Integer.toString(left + ddx));
            border_2.setAttributeNS(null, "y2", Integer.toString(gridSize + ddy));
            border_2.setAttributeNS(null, "stroke", "gray");
            border_2.setAttributeNS(null, "stroke-width", Integer.toString(borderLineSize));
            svgRoot.appendChild(border_2);

            Element border_3 = doc.createElementNS(svgNS, "line");
            border_3.setAttributeNS(null, "x1", Integer.toString(ddx));
            border_3.setAttributeNS(null, "y1", Integer.toString(right + ddy));
            border_3.setAttributeNS(null, "x2", Integer.toString(gridSize + ddx));
            border_3.setAttributeNS(null, "y2", Integer.toString(right + ddy));
            border_3.setAttributeNS(null, "stroke", "gray");
            border_3.setAttributeNS(null, "stroke-width", Integer.toString(borderLineSize));
            svgRoot.appendChild(border_3);

            Element border_4 = doc.createElementNS(svgNS, "line");
            border_4.setAttributeNS(null, "x1", Integer.toString(right + ddx));
            border_4.setAttributeNS(null, "y1", Integer.toString(ddy));
            border_4.setAttributeNS(null, "x2", Integer.toString(right + ddx));
            border_4.setAttributeNS(null, "y2", Integer.toString(gridSize + ddy));
            border_4.setAttributeNS(null, "stroke", "gray");
            border_4.setAttributeNS(null, "stroke-width", Integer.toString(borderLineSize));
            svgRoot.appendChild(border_4);
        }

        if (drawInnerLines) { // grid
            for (int i = 1; i < n; i++) {
                int coord = i * (cellSize + borderLineSize) + borderLineSize / 2;

                Element inner_line_1 = doc.createElementNS(svgNS, "line");
                inner_line_1.setAttributeNS(null, "x1", Integer.toString(borderLineSize + ddx));
                inner_line_1.setAttributeNS(null, "y1", Integer.toString(coord + ddy));
                inner_line_1.setAttributeNS(null, "x2", Integer.toString(gridSize - borderLineSize + ddx));
                inner_line_1.setAttributeNS(null, "y2", Integer.toString(coord + ddy));
                inner_line_1.setAttributeNS(null, "stroke", "gray");
                inner_line_1.setAttributeNS(null, "stroke-width", Integer.toString(borderLineSize));
                svgRoot.appendChild(inner_line_1);

                Element inner_line_2 = doc.createElementNS(svgNS, "line");
                inner_line_2.setAttributeNS(null, "x1", Integer.toString(coord + ddx));
                inner_line_2.setAttributeNS(null, "y1", Integer.toString(borderLineSize + ddy));
                inner_line_2.setAttributeNS(null, "x2", Integer.toString(coord + ddx));
                inner_line_2.setAttributeNS(null, "y2", Integer.toString(gridSize - borderLineSize + ddy));
                inner_line_2.setAttributeNS(null, "stroke", "gray");
                inner_line_2.setAttributeNS(null, "stroke-width", Integer.toString(borderLineSize));
                svgRoot.appendChild(inner_line_2);
            }
        }

        // filling
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Element rectangle = doc.createElementNS(svgNS, "rect");
                rectangle.setAttributeNS(null, "x", 
                        Integer.toString(ddx + i * cellSize + (i + 1) * borderLineSize));
                rectangle.setAttributeNS(null, "y",
                        Integer.toString(ddy + j * cellSize + (j + 1) * borderLineSize));
                rectangle.setAttributeNS(null, "width", Integer.toString(cellSize));
                rectangle.setAttributeNS(null, "height", Integer.toString(cellSize));
                rectangle.setAttributeNS(null, "fill", getColorXML(distMatrix[perm[j]][perm[i]]));
                svgRoot.appendChild(rectangle);
            }
        }

        return doc;
    }

    public Document createHeatMapWithNames(Document doc, int ddx, int ddy) {
        int width = gridSize + dx_for_names; // +dx for row names
        int height = gridSize + dy_for_names;

        Element svgRoot = doc.getDocumentElement();

        // drawing heat map itself
        Document heatMap = createHeatMap(ddx, ddy + dy_for_names);
        NodeList nodes = heatMap.getDocumentElement().getChildNodes();
        for (int i = 0; i < nodes.getLength(); i++) {
            org.w3c.dom.Node node = doc.importNode(nodes.item(i), false);
            svgRoot.appendChild(node);
        }

        double theta = -90;
        double anchorx = ddx + gridSize / 2.0;
        double anchory = ddy + dy_for_names + gridSize / 2.0;
        String transform = "translate(" + anchorx + "," + anchory + 
                            ") rotate(" + theta + ") translate( " + -anchorx + "," + -anchory + ")";
        // drawing row and column names
        for (int i = 0; i < n; i++) {
            int yc = dy_for_names + (i + 1) * (cellSize + borderLineSize) - cellSize / 3;
            Element rowName = doc.createElementNS(svgNS, "text");
            rowName.setAttributeNS(null, "x", Integer.toString(ddx + gridSize + 10));
            rowName.setAttributeNS(null, "y", Integer.toString(ddy + yc));
            rowName.setAttributeNS(null, "font", fontXML);
            rowName.setAttributeNS(null, "fill", "black");
            rowName.setTextContent(names[perm[i]]);
            svgRoot.appendChild(rowName);

            Element columnName = doc.createElementNS(svgNS, "text");
            columnName.setAttributeNS(null, "x", Integer.toString(ddx + gridSize + 10));
            columnName.setAttributeNS(null, "y", Integer.toString(ddy + yc));
            columnName.setAttributeNS(null, "font", fontXML);
            columnName.setAttributeNS(null, "fill", "black");
            columnName.setAttributeNS(null, "transform", transform);
            columnName.setTextContent(names[perm[i]]);
            svgRoot.appendChild(columnName);
        }

        return doc;
    }

    private String getColorXML(double value) {
        double p = (value - low) / (high - low);
        int r = Math.min(255,
                lowColor.getRed() + (int) Math.round(p * (highColor.getRed() - lowColor.getRed())));
        int g = Math.min(255,
                lowColor.getGreen() + (int) Math.round(p * (highColor.getGreen() - lowColor.getGreen())));
        int b = Math.min(255,
                lowColor.getBlue() + (int) Math.round(p * (highColor.getBlue() - lowColor.getBlue())));
        return "rgb(" + r + "," + g + "," + b + ")";
    }


    // ============================   Clustering objects   ============================

    class Node {
        int no = -1;    // object number, or -1 - if not a leaf
        int leafs = 1;  // number of leafs

        Node left, right;
        double distance;

        int dy, dx; // dy - from top to bottom, dx - from right to left
    }


    /**
     * Clusters objects, works in O(n^3) time.
     * @return the resulting root node.
     */
    public Node clusterObjects() {
        // preparing
        Node[] nodes = new Node[n];
        for (int i = 0; i < n; i++) {
            Node c = new Node();
            c.no = i;
            c.dx = 0;
            c.dy = (i + 1) * borderLineSize + i * cellSize + cellSize / 2;
            nodes[i] = c;
        }

        double[][] dist = new double[n][n];
        int[] g1 = new int[n];
        int[] g2 = new int[n];

        for (int i = 0; i < n; i++) {
            getGroup(nodes[i], g1, 0);
            dist[i][i] = 0;
            for (int j = i + 1; j < n; j++) {
                getGroup(nodes[j], g2, 0);
                dist[i][j] = dist[j][i] = distanceBetweenGroups(g1, 1, g2, 1);
            }
        }

        // clustering
        int count = n;
        Node root = (n > 0) ? nodes[0] : null;
        while (count > 1) {
//            System.out.println("Count = " + count);

            double minDist = Double.MAX_VALUE;
            int i = -1, j = -1;

            for (int ii = 0; ii < n; ii++) {
                for (int jj = ii + 1; jj < n; jj++) {
                    if (nodes[ii] != null && nodes[jj] != null && dist[ii][jj] < minDist) {
                        minDist = dist[ii][jj];
                        i = ii;
                        j = jj;
                    }
                }
            }

            if (i == -1 || minDist < 0) {
                throw new RuntimeException("Internal error. Wrong minDist index.");
            }

            // merging node i and j
//            System.out.println("Merging " + i + " and " + j + ", min dist = " + minDist);
            root = new Node();
            root.left = nodes[i];
            root.right = nodes[j];
            root.distance = minDist;
            root.leafs = root.left.leafs + root.right.leafs;

            root.dy = (root.left.dy + root.right.dy) / 2;
            root.dx = (int) (minDist / high * dx_for_dendrogram);

            // updating
            nodes[i] = root;
            nodes[j] = null;
            int g1n = getGroup(root, g1, 0);
            for (int ii = 0; ii < n; ii++) {
                dist[ii][j] = dist[j][ii] = -1;
                if (ii != i) {
                    int g2n = getGroup(nodes[ii], g2, 0);
                    dist[ii][i] = dist[i][ii] = distanceBetweenGroups(g1, g1n, g2, g2n);
                }
            }
            count--;
        }

//        System.out.println("Done");
        return root;
    }

    protected double distanceBetweenGroups(int[] g1, int g1n, int[] g2, int g2n) {
        if (g1n == 0 || g2n == 0) {
            return -1;
        }
        double sum = 0;
        for (int i = 0; i < g1n; i++) {
            for (int j = 0; j < g2n; j++) {
                sum += distMatrix[g1[i]][g2[j]];
            }
        }
        return sum / g1n / g2n;
    }

    int getGroup(Node n, int[] e, int first) {
        if (n == null) {
            return 0;
        }
        if (n.no >= 0) {
            e[first] = n.no;
            return 1;
        }
        int c = getGroup(n.left,  e, first);
        c    += getGroup(n.right, e, first + c);
        return c;
    }

    protected void renumber(Node node, int first) {
        if (node.no >= 0) {
            // node is leaf
            perm[first] = node.no;
            node.dy = (first + 1) * borderLineSize + first * cellSize + cellSize / 2;
            return;
        }

        renumber(node.left, first);
        renumber(node.right, first + node.left.leafs);
        node.dy = (node.left.dy + node.right.dy) / 2;
    }


    public Document createLeftDendrogram(Document doc, boolean shouldRenumber, int ddx, int ddy) {
        int width = dx_before_dendrogram + dx_for_dendrogram;
        int height = gridSize;

        Element svgRoot = doc.getDocumentElement();

        Node root = clusterObjects();
        if (shouldRenumber) {
            renumber(root, 0);
        }

        drawClusterNode(root, doc, width, ddx, ddy);
        
        return doc;
    }

    void drawClusterNode(Node n, Document doc, int width, int ddx, int ddy) {
        if (n == null || n.no >= 0) {
            return;
        }

        Element line_1 = doc.createElementNS(svgNS, "line");
        line_1.setAttributeNS(null, "x1", Integer.toString(ddx + width - n.left.dx));
        line_1.setAttributeNS(null, "y1", Integer.toString(ddy + n.left.dy));
        line_1.setAttributeNS(null, "x2", Integer.toString(ddx + width - n.dx));
        line_1.setAttributeNS(null, "y2", Integer.toString(ddy + n.left.dy));
        line_1.setAttributeNS(null, "stroke", "black");
        line_1.setAttributeNS(null, "stroke-width", Integer.toString(borderLineSize));
        doc.getDocumentElement().appendChild(line_1);

        Element line_2 = doc.createElementNS(svgNS, "line");
        line_2.setAttributeNS(null, "x1", Integer.toString(ddx + width - n.right.dx));
        line_2.setAttributeNS(null, "y1", Integer.toString(ddy + n.right.dy));
        line_2.setAttributeNS(null, "x2", Integer.toString(ddx + width - n.dx));
        line_2.setAttributeNS(null, "y2", Integer.toString(ddy + n.right.dy));
        line_2.setAttributeNS(null, "stroke", "black");
        line_2.setAttributeNS(null, "stroke-width", Integer.toString(borderLineSize));
        doc.getDocumentElement().appendChild(line_2);

        Element line_3 = doc.createElementNS(svgNS, "line");
        line_3.setAttributeNS(null, "x1", Integer.toString(ddx + width - n.dx));
        line_3.setAttributeNS(null, "y1", Integer.toString(ddy + n.left.dy));
        line_3.setAttributeNS(null, "x2", Integer.toString(ddx + width - n.dx));
        line_3.setAttributeNS(null, "y2", Integer.toString(ddy + n.right.dy));
        line_3.setAttributeNS(null, "stroke", "black");
        line_3.setAttributeNS(null, "stroke-width", Integer.toString(borderLineSize));
        doc.getDocumentElement().appendChild(line_3);

        drawClusterNode(n.left, doc, width, ddx, ddy);
        drawClusterNode(n.right, doc, width, ddx, ddy);

    }



    // ===========================   Drawing full heat map   ================================

    public Document createColorScale(Document doc, int ddx, int ddy) {
        int dy_before = 20;
        int dy_after = 30;

        Element svgRoot = doc.getDocumentElement();

        int dy_dist = dy_before + dy_scale + 20;
        { // color distance
            Element text_1 = doc.createElementNS(svgNS, "text");
            text_1.setAttributeNS(null, "x", Integer.toString(ddx + 10));
            text_1.setAttributeNS(null, "y", Integer.toString(ddy + dy_before + dy_scale / 2 + 10));
            text_1.setAttributeNS(null, "font", fontXML);
            text_1.setAttributeNS(null, "fill", "black");
            text_1.setTextContent("Color");
            svgRoot.appendChild(text_1);

            Element text_2 = doc.createElementNS(svgNS, "text");
            text_2.setAttributeNS(null, "x", Integer.toString(ddx + 10));
            text_2.setAttributeNS(null, "y", Integer.toString(ddy + dy_dist));
            text_2.setAttributeNS(null, "font", fontXML);
            text_2.setAttributeNS(null, "fill", "black");
            text_2.setTextContent("Distance");
            svgRoot.appendChild(text_2);
        }

        // colors and labels
        int cellSize = dx_scale / 6;
        for (int i = 0; i <= 5; i++) {
            double v = low + (high - low) * i / 5.0;

            Element rectangle = doc.createElementNS(svgNS, "rect");
            rectangle.setAttributeNS(null, "x", Integer.toString(ddx + dx_scale_text + i * cellSize));
            rectangle.setAttributeNS(null, "y", Integer.toString(ddy + dy_before));
            rectangle.setAttributeNS(null, "width", Integer.toString(cellSize));
            rectangle.setAttributeNS(null, "height", Integer.toString(dy_scale));
            rectangle.setAttributeNS(null, "fill", getColorXML(v));
            svgRoot.appendChild(rectangle);

            Element text = doc.createElementNS(svgNS, "text");
            text.setAttributeNS(null, "x", Integer.toString(ddx + dx_scale_text + i * cellSize + cellSize / 5));
            text.setAttributeNS(null, "y", Integer.toString(ddy + dy_dist));
            text.setAttributeNS(null, "font", fontXML);
            text.setAttributeNS(null, "fill", "black");
            text.setTextContent(String.format("%.1f", v));
            svgRoot.appendChild(text);
        }

        return doc;
    }

    public Document createFullHeatMap(boolean shouldRenumber) {
        int width = dx_before_dendrogram + dx_for_dendrogram + gridSize + dx_for_names;
        int height = dy_for_names + gridSize + dy_before_color + dy_scale + dy_after_color;

        DOMImplementation impl = SVGDOMImplementation.getDOMImplementation();
        Document doc = impl.createDocument(svgNS, "svg", null);

        Element svgRoot = doc.getDocumentElement();
        svgRoot.setAttributeNS(null, "width", Integer.toString(width));
        svgRoot.setAttributeNS(null, "height", Integer.toString(height));

        { // background
            Element rectangle = doc.createElementNS(svgNS, "rect");
            rectangle.setAttributeNS(null, "x", Integer.toString(0));
            rectangle.setAttributeNS(null, "y", Integer.toString(0));
            rectangle.setAttributeNS(null, "width", Integer.toString(width));
            rectangle.setAttributeNS(null, "height", Integer.toString(height));
            rectangle.setAttributeNS(null, "fill", "white");
            svgRoot.appendChild(rectangle);
        }

        if (n >= 2) {
            createLeftDendrogram(doc, shouldRenumber, 0, dy_for_names);
        }

        createHeatMapWithNames(doc, dx_before_dendrogram + dx_for_dendrogram, 0);

        createColorScale(doc, dx_before_dendrogram + dx_for_dendrogram - 120,
                dy_for_names + gridSize + dy_before_color);

        return doc;
    }


    public static void main(String[] args) throws IOException {
        double[][] matrix = {{  0, 0.9,   1, 0.9, 0.8},
                             {0.9,   0, 0.5, 0.6, 0.1},
                             {  1, 0.5,   0, 0.2, 0.6},
                             {0.9, 0.6, 0.2,   0, 0.7},
                             {0.8, 0.1, 0.6, 0.7, 0}};
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                for (int k = 0; k < matrix.length; k++) {
//                    if (i != j && i != k && j != k) {
                        if (matrix[i][j] + matrix[j][k] < matrix[i][k]) {
                            System.out.println("Wrong for i = " + i + ", j = " + j + ", k = " + k);
                            System.out.println((matrix[i][j] + matrix[j][k]) + " < " + matrix[i][k]);
                            return;
                        }
//                    }
                }
            }
        }
        String[] names = {"Abracadabra library got from some url", "Boldii fish lib2s",
                "Cnot", "Don't know", "mmmm"};

        FullHeatMapXML hm = new FullHeatMapXML(matrix, names);
        OutputFormat format = new OutputFormat(Method.XML, "utf-8", true);
        format.setIndent(4);
        format.setLineWidth(80);
        format.setPreserveEmptyAttributes(true);
        format.setPreserveSpace(true);
        XMLSerializer xmlSerializer = new XMLSerializer(new FileWriter("test.svg"), format);
        xmlSerializer.serialize(hm.createFullHeatMap(true));
    }

}
