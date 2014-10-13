public class ToolsScanner extends ru.ifmo.genetics.ToolsScanner {

    public ToolsScanner(String toolsPath, Class classInYourPackage) {
        super(toolsPath, classInYourPackage);
    }


    public static void main(String[] args) {
        new ToolsScanner(
                "tools", ToolsScanner.class
        ).run();
    }

}
