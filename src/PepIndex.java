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

    private ArrayList<Path> peptideFiles;

    @Override
    public FileVisitResult preVisitDirectory(Path dir,
        BasicFileAttributes attrs) throws IOException {
      System.out.println("Processing directory " + dir);
      peptideFiles = new ArrayList<>();
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

      if (textFileMatcher.matches(file.getFileName())) peptideFiles.add(file);
      else if (fastaFileMatcher.matches(file.getFileName())) {
        readProteins(file);
      } else {
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
        System.out.println("No protein fasta file found so far." + dir);
      else for (Path peptideFile : peptideFiles)
        findPeptidePositions(peptideFile);

      return FileVisitResult.CONTINUE;
    }

  }

  private static void readProteins(Path proteinFile) {
    proteinIdToSeq = new HashMap<>();
    try (Scanner scan = new Scanner(proteinFile)) {
      String line = null;
      if (scan.hasNextLine()) line = scan.nextLine();
      while (line != null && line.startsWith(">")) {
        int proteinIdStartIndex = line.indexOf("|");
        if (proteinIdStartIndex == -1)
          System.out.println("Cannot find | in line starting with >");
        String proteinID = line.substring(proteinIdStartIndex + 1,
            line.indexOf("|", proteinIdStartIndex + 1));
        StringBuilder proteinPeptides = new StringBuilder();
        while ((line = scan.nextLine()) != null && line.length() > 0
            && !line.startsWith(">"))
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
    final String delimiter = "\t";
    StringBuilder columns = new StringBuilder();

    int nameCount = peptideFile.getNameCount();

    StringBuilder currentOutputDir = new StringBuilder("output/");

    if (nameCount > 2) {

      for (int i = 1; i < nameCount - 1; i++)
        currentOutputDir.append(peptideFile.getName(i));

      try {
        Files.createDirectories(Paths.get(currentOutputDir.toString()));
      } catch (IOException e2) {
        System.out.println("Cannot properly create dir for " + peptideFile);
        System.exit(-1);
      }
    }

    try (BufferedReader br = Files.newBufferedReader(peptideFile)) {
      try (
          BufferedWriter bw = Files.newBufferedWriter(
              Paths.get(currentOutputDir.toString() + "/"
                  + peptideFile.getFileName()),
              StandardCharsets.UTF_8, StandardOpenOption.CREATE,
              StandardOpenOption.WRITE)) {

        int count = 0;

        String headerLine = br.readLine();

        String line;
        while ((line = br.readLine()) != null) {

          ArrayList<String> peptides = new ArrayList<>();
          ArrayList<Integer> peptideShifts = new ArrayList<>();
          ArrayList<String> proteinIds = new ArrayList<>();

          int index = line.indexOf("-.");
          peptides.add(line.substring(index + 2,
              (index = line.indexOf("(", index + 2))));
          peptideShifts.add(Integer
              .valueOf(line.substring(index + 1, line.indexOf(")", index))));

          index = line.indexOf("--");
          peptides.add(line.substring(index + 2,
              (index = line.indexOf("(", index + 2))));
          peptideShifts.add(Integer
              .valueOf(line.substring(index + 1, line.indexOf(")", index))));

          while ((index = line.indexOf("sp|", index)) != -1) {
            proteinIds.add(line.substring(index + 3,
                (index = line.indexOf("|", index + 3))));
          }

          for (int i = 1; i < proteinIds.size(); i++)
            columns.append(delimiter).append("protein").append(i)
                .append("info");

          for (int i = 1; i <= peptides.size(); i++)
            columns.append(delimiter).append("peptide").append(i)
                .append(delimiter).append("position").append(i)
                .append(delimiter).append("proteinId").append(i);

          if (count == 0) {
            bw.write(headerLine + columns);
            bw.newLine();
          }

          bw.write(line);

          for (int i = 0; i < peptides.size(); i++) {
            String peptide = peptides.get(i);
            ArrayList<Integer> result = findPosition(peptide, proteinIds);
            bw.write(delimiter + peptide + delimiter
                + (result.get(1) + peptideShifts.get(i)));
            bw.write(delimiter + proteinIds.get(result.get(0)));
          }

          bw.newLine();

          count++;
        }

      } catch (IOException e) {
        System.out.println("Cannot write output file: " + currentOutputDir + "/"
            + peptideFile.getFileName());
        System.exit(-1);
      }
    } catch (IOException e1) {
      System.out.println("Cannot read peptide file " + peptideFile);
      System.exit(-1);
    }
  }

  private static ArrayList<Integer> findPosition(String peptide,
      ArrayList<String> proteinIds) {

    int position = -1;
    ArrayList<Integer> result = new ArrayList<>();

    for (int i = 0; i < proteinIds.size(); i++) {
      int temp = proteinIdToSeq.get(proteinIds.get(i)).indexOf(peptide);
      // System.out.println("protein " + protein + " peptide " + peptide);
      if (position != -1 && temp != -1) {
        System.out.printf("Find peptide %s in both proteins.%n", peptide);
        System.exit(-1);
      }
      if (temp != -1) {
        position = temp;
        result.add(i);
      }
    }

    if (position == -1) {
      System.out.println("Cannot find peptide " + peptide);
      System.exit(-1);
    }

    result.add(position);

    return result;
  }

  public static void main(String[] args) {

    try {
      Files.walkFileTree(Paths.get(INPUT_DIR), new ProcessInputFiles());
    } catch (IOException e) {
      System.out.println("Cannot read input directory: " + INPUT_DIR);
    }

  }
}
