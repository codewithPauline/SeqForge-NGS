import java.nio.file.Files
import java.nio.file.Path

class Inputs {
    static Path localFile(Object value, Path base, String label) {
        if (!value || value.toString().contains('://'))
            throw new IllegalArgumentException("${label}: provide an existing local file path")
        Path path = base.resolve(value.toString().trim()).normalize().toAbsolutePath()
        if (!Files.isRegularFile(path) || !Files.isReadable(path) || Files.size(path) == 0)
            throw new IllegalArgumentException("${label}: file is missing, empty, or unreadable: ${path}")
        return path.toRealPath()
    }

    static List samples(List rows, Path base) {
        if (!rows) throw new IllegalArgumentException('Samplesheet contains no samples')
        def ids = [] as Set
        def reads = [] as Set
        return rows.collect { row ->
            if (row.keySet() != ['sample', 'fastq_1', 'fastq_2'] as Set)
                throw new IllegalArgumentException('Samplesheet must have exactly: sample,fastq_1,fastq_2')
            def id = row.sample?.toString()?.trim()
            if (!id || !(id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]{0,79}/))
                throw new IllegalArgumentException("Invalid sample ID: ${id}; use 1–80 letters, digits, dots, underscores, or hyphens")
            if (!ids.add(id)) throw new IllegalArgumentException("Duplicate sample ID: ${id}; one read pair per sample is supported")
            def mates = ['fastq_1', 'fastq_2'].collect { column ->
                def path = localFile(row[column], base, "${id} ${column}")
                if (!(path.fileName.toString() ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/))
                    throw new IllegalArgumentException("${id}: FASTQ filenames must use letters, digits, dots, underscores, or hyphens")
                if (!(path.toString() ==~ /.*\.(fastq|fq)(\.gz)?/))
                    throw new IllegalArgumentException("${id}: expected .fastq/.fq, optionally gzipped: ${path}")
                if (!reads.add(path)) throw new IllegalArgumentException("FASTQ reused in samplesheet: ${path}")
                return path
            }
            return [id, mates[0], mates[1]]
        }
    }
}
