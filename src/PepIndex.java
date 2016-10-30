import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.FileSystems;
import java.nio.file.FileVisitResult;
import java.nio.file.FileVisitor;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;
import java.util.Scanner;

/**
 * Find indexes for peptide 1 and peptide 2.
 * 
 * @author Zexi Jesse Zhuang
 */
public class PepIndex {
  /**
   * map protein Id to protein peptide sequence.
   */
  private static Map<String, String> proteinIdToSeq;

  private final static String INPUT_DIR = "input/";

  private static class ProcessInputFiles implements FileVisitor<Path> {

    private ArrayList<Path> dataFiles;

    public ProcessInputFiles() {
      dataFiles = new ArrayList<>();
    }

    @Override
    public FileVisitResult preVisitDirectory(Path dir,
        BasicFileAttributes attrs) throws IOException {
      System.out.println("Processing directory " + dir);
      Objects.requireNonNull(dir);
      Objects.requireNonNull(attrs);
      return FileVisitResult.CONTINUE;
    }

    @Override
    public FileVisitResult visitFile(Path file, BasicFileAttributes attrs)
        throws IOException {
      Objects.requireNonNull(file);
      Objects.requireNonNull(attrs);
      System.out.println(">>>Processing file " + file);

      PathMatcher textFileMatcher = FileSystems.getDefault()
          .getPathMatcher("glob:*.txt");
      PathMatcher fastaFileMatcher = FileSystems.getDefault()
          .getPathMatcher("glob:*.fasta");

      if (textFileMatcher.matches(file.getFileName())) dataFiles.add(file);
      else if (fastaFileMatcher.matches(file.getFileName())) readProteins(file);
      else {
        System.out.println("Unknown file type: " + file);
        System.exit(-1);
      }

      return FileVisitResult.CONTINUE;
    }

    @Override
    public FileVisitResult visitFileFailed(Path file, IOException exc)
        throws IOException {
      Objects.requireNonNull(file);
      System.out.println("Cannot visit file " + file);
      throw exc;
    }

    @Override
    public FileVisitResult postVisitDirectory(Path dir, IOException exc)
        throws IOException {
      System.out.println("Finished processing directory " + dir);
      Objects.requireNonNull(dir);
      if (exc != null) throw exc;

      if (proteinIdToSeq == null)
        System.out.println("No protein fasta file found.");
      else for (Path peptideFile : dataFiles)
        findPeptidePositions(peptideFile);

      return FileVisitResult.CONTINUE;
    }

  }

  private static void readProteins(Path proteinFile) {
    proteinIdToSeq = new HashMap<>();
    try (Scanner scan = new Scanner(proteinFile)) {
      String firstLine;
      while (scan.hasNextLine() && (firstLine = scan.nextLine()) != null
          && firstLine.startsWith(">")) {
        int proteinIdStartIndex = firstLine.indexOf("|");
        if (proteinIdStartIndex == -1)
          System.out.println("Cannot find | in line starting with >");
        String proteinID = firstLine.substring(proteinIdStartIndex + 1,
            firstLine.indexOf("|", proteinIdStartIndex + 1));
        StringBuilder proteinPeptides = new StringBuilder();
        String line;
        while ((line = scan.nextLine()) != null && line.length() > 0)
          proteinPeptides.append(line);
        proteinIdToSeq.put(proteinID, proteinPeptides.toString());
      }

    } catch (IOException e) {
      System.out.println("Cannot open protein file: " + proteinFile);
      System.exit(-1);
    }
    System.out.println(proteinIdToSeq);
  }

  private static void findPeptidePositions(Path peptideFile) {
    final String outputDirectory = "output/";
    final String delimiter = "\t";
    final String columns = delimiter + "peptide1" + delimiter + "position1"
        + delimiter + "peptide2" + delimiter + "position2";

    try (BufferedReader br = Files.newBufferedReader(peptideFile)) {
      try (BufferedWriter bw = Files.newBufferedWriter(
          Paths.get(outputDirectory + peptideFile.getFileName()),
          StandardCharsets.UTF_8, StandardOpenOption.CREATE,
          StandardOpenOption.WRITE)) {

        bw.write(br.readLine() + columns);
        bw.newLine();

        String line;
        while ((line = br.readLine()) != null) {
          int index = line.indexOf("-.");
          String peptide1 = line.substring(index + 2,
              (index = line.indexOf("(", index + 2)));

          int peptide1Pos = Integer
              .valueOf(line.substring(index + 1, line.indexOf(")", index)));
          index = line.indexOf("--");
          String peptide2 = line.substring(index + 2,
              (index = line.indexOf("(", index + 2)));
          int peptide2Pos = Integer
              .valueOf(line.substring(index + 1, line.indexOf(")", index)));

          index = line.indexOf("sp|");
          final String proteinId = line.substring(index + 3,
              line.indexOf("|", index + 3));
          final String proteinPeptides = proteinIdToSeq.get(proteinId);

          bw.write(line + delimiter + peptide1 + delimiter
              + (proteinPeptides.indexOf(peptide1) + peptide1Pos) + delimiter
              + peptide2 + delimiter
              + (proteinPeptides.indexOf(peptide2) + peptide2Pos));
          bw.newLine();
        }

      } catch (IOException e) {
        System.out.println(
            "Cannot write output file: " + outputDirectory + peptideFile);
        System.exit(-1);
      }
    } catch (IOException e1) {
      System.out.println("Cannot read peptide file " + INPUT_DIR + peptideFile);
      System.exit(-1);
    }
  }

  public static void main(String[] args) {

    // Path path = Paths
    // .get("input/20161021_ProAI_8nm_EDC_01.validated.intra.txt");
    // PathMatcher matcher =
    // FileSystems.getDefault().getPathMatcher("glob:*.txt");
    // System.out.println(path.getFileName());
    // System.out.println(path.getNameCount());
    // System.out.println(path.getName(1));
    // System.out.println(matcher.matches(path.getFileName()));

    try {
      Files.walkFileTree(Paths.get(INPUT_DIR), new ProcessInputFiles());
    } catch (IOException e) {
      System.out.println("Cannot read input directory: " + INPUT_DIR + "/");
    }

  }
}
