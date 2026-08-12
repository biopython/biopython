class AlignmentIterator(bigbed.AlignmentIterator, maf.AlignmentIterator):
    """Alignment iterator for bigMaf files.

    The file may contain multiple alignments, which are loaded and returned
    incrementally.

    Alignment annotations are stored in the ``.annotations`` attribute of
    the ``Alignment`` object, except for the alignment score, which is stored
    as an attribute. Sequence information of empty parts in the alignment
    block is stored in the alignment annotations under the ``"empty"`` key.
    Annotations specific to each line in the alignment are stored in the
    ``.annotations`` attribute of the corresponding sequence record.
    """

    fmt = "bigMaf"
    mode = "b"

    def __init__(self, source):
        """Create an AlignmentIterator object."""
        self.reference = None
        self._index = 0
        super().__init__(source)

    def _read_reference(self, stream):
        """Read the reference sequence name from the first alignment."""
        formatter = struct.Struct(self.byteorder + "III")
        size = formatter.size

        node = self.tree
        while True:
            try:
                node = node.children[0]
            except AttributeError:
                break

        file_position = stream.tell()
        stream.seek(node.dataOffset)

        data_size = 256
        data = b""
        compressed_data = b""

        while True:
            chunk = stream.read(data_size)

            if self._compressed:
                compressed_data += chunk
                decompressor = zlib.decompressobj()
                data = decompressor.decompress(compressed_data)
            else:
                data += chunk

            try:
                index = data.index(b";s", size)
            except ValueError:
                continue

            words = data[index + 1:].split()
            if len(words) > 2:
                break

        name = words[1]
        stream.seek(file_position)

        reference, _chromosome = name.split(b".", 1)

        return reference.decode()

    def _read_header(self, stream):
        """Read the bigMaf header."""
        super()._read_header(stream)

        if self.reference is None:
            self.reference = self._read_reference(stream)

        self.targets[0].id = f"{self.reference}.{self.targets[0].id}"

    @staticmethod
    def _read_line(buffer, index):
        """Read a complete line from the alignment buffer."""
        match = re.match(b"^[^;]*", buffer[index:])
        return index + match.span()[1]

    def _parse_annotation(self, buffer, index, annotations):
        """Parse an MAF annotation line."""
        end = self._read_line(buffer, index)
        line = buffer[index:end].tobytes()
        words = line[1:].split()

        score = None

        for word in words:
            key, value = word.split(b"=")

            if key == b"score":
                score = float(value)

            elif key == b"pass":
                pass_value = int(value)

                if pass_value <= 0:
                    raise ValueError(
                        f"pass value must be positive (found {pass_value})"
                    )

                annotations["pass"] = pass_value

            else:
                raise ValueError(
                    f"Unknown annotation variable '{key.decode()}'"
                )

        return end, score

    def _parse_sequence(
        self,
        buffer,
        index,
        data,
        data_start,
        printed_alignment_parser,
        records,
        strands,
    ):
        """Parse an MAF sequence (s) line."""
        match = re.match(
            rb"^s\s*\S*\s*\d*\s*\d*\s*[+-]\s*\d*\s*",
            buffer[index:],
        )
        end = index + match.span()[1]

        line = buffer[index:end].tobytes()
        words = line.split(None, 5)

        if len(words) != 6:
            raise ValueError(
                "Error parsing alignment - 's' line must have 7 fields"
            )

        src = words[1].decode()
        start = int(words[2])
        sequence_size = int(words[3])
        strand = words[4]
        src_size = int(words[5])

        sequence_start = end

        n, sequence = printed_alignment_parser.feed(
            data,
            data_start + sequence_start,
        )

        if len(sequence) != sequence_size:
            raise ValueError(
                "sequence size is incorrect "
                f"(found {len(sequence)}, expected {sequence_size})"
            )

        seq = Seq({start: sequence}, length=src_size)

        record = SeqRecord(
            seq,
            id=src,
            name="",
            description="",
        )

        records.append(record)
        strands.append(strand)

        return sequence_start + n, src, record

    def _parse_insert(self, buffer, index, src, record):
        """Parse an MAF insertion (i) line."""
        end = self._read_line(buffer, index)
        line = buffer[index:end].tobytes()
        words = line.split(None, 5)

        assert len(words) == 6
        assert words[1].decode() == src

        left_status = words[2].decode()
        left_count = int(words[3])
        right_status = words[4].decode()
        right_count = int(words[5])

        assert left_status in self.status_characters
        assert right_status in self.status_characters

        record.annotations["leftStatus"] = left_status
        record.annotations["leftCount"] = left_count
        record.annotations["rightStatus"] = right_status
        record.annotations["rightCount"] = right_count

        return end

    def _parse_empty(self, buffer, index, annotations):
        """Parse an MAF empty (e) line."""
        end = self._read_line(buffer, index)
        line = buffer[index:end].tobytes()
        words = line.split(None, 6)

        assert len(words) == 7

        src = words[1].decode()
        start = int(words[2])
        sequence_size = int(words[3])
        strand = words[4]
        src_size = int(words[5])
        status = words[6].decode()

        assert status in self.empty_status_characters

        sequence = Seq(None, length=src_size)

        record = SeqRecord(
            sequence,
            id=src,
            name="",
            description="",
        )

        end_position = start + sequence_size

        if strand == b"+":
            segment = (start, end_position)
        else:
            segment = (
                src_size - start,
                src_size - end_position,
            )

        empty = (record, segment, status)

        annotations.setdefault("empty", []).append(empty)

        return end, src

    def _parse_quality(self, buffer, index, src, record):
        """Parse an MAF quality (q) line."""
        end = self._read_line(buffer, index)
        line = buffer[index:end].tobytes()
        words = line.split(None, 2)

        assert len(words) == 3
        assert words[1].decode() == src

        value = words[2].replace(b"-", b"")
        record.annotations["quality"] = value.decode()

        return end

    @staticmethod
    def _finish_alignment(
        records,
        strands,
        annotations,
        score,
        printed_alignment_parser,
    ):
        """Build the final Alignment object."""
        coordinates = np.empty(
            printed_alignment_parser.shape,
            np.intp,
        )

        printed_alignment_parser.fill(coordinates)

        for record, strand, row in zip(
            records,
            strands,
            coordinates,
        ):
            if strand == b"+":
                row += record.seq.defined_ranges[0][0]
            else:
                record.seq = record.seq.reverse_complement()
                row[:] = record.seq.defined_ranges[0][1] - row

        alignment = Alignment(records, coordinates)

        if annotations:
            alignment.annotations = annotations

        if score is not None:
            alignment.score = score

        return alignment

    def _parse_block(
        self,
        buffer,
        data,
        data_start,
    ):
        """Parse all lines from one MAF block."""
        records = []
        strands = []
        annotations = {}
        score = None

        parser = _aligncore.PrintedAlignmentParser(b";")

        index = -1
        src = None
        record = None

        while True:
            index += 1
            prefix = buffer[index:index + 1]

            if prefix == b"#":
                index = self._read_line(buffer, index) - 1

            elif prefix == b"a":
                index, parsed_score = self._parse_annotation(
                    buffer,
                    index,
                    annotations,
                )

                if parsed_score is not None:
                    score = parsed_score

                index -= 1

            elif prefix == b"s":
                index, src, record = self._parse_sequence(
                    buffer,
                    index,
                    data,
                    data_start,
                    parser,
                    records,
                    strands,
                )

                index -= 1

            elif prefix == b"i":
                index = self._parse_insert(
                    buffer,
                    index,
                    src,
                    record,
                ) - 1

            elif prefix == b"e":
                index, src = self._parse_empty(
                    buffer,
                    index,
                    annotations,
                )

                index -= 1

            elif prefix == b"q":
                index = self._parse_quality(
                    buffer,
                    index,
                    src,
                    record,
                ) - 1

            elif prefix == b"\00":
                break

            else:
                line = buffer[index:].tobytes()
                raise ValueError(
                    "Error parsing alignment - "
                    f"unexpected line:\n{line}"
                )

        return (
            records,
            strands,
            annotations,
            score,
            parser,
        )

    # pylint: disable=invalid-name, too-many-arguments
    def _create_alignment(
        self,
        chromId,
        chromStart,
        chromEnd,
        data,
        dataStart,
        dataEnd,
    ):
        """Create an Alignment object from a bigMaf data block."""
        del chromId, chromStart, chromEnd

        assert data[dataEnd - 1] == 0

        buffer = memoryview(data)[dataStart:dataEnd]

        (
            records,
            strands,
            annotations,
            score,
            parser,
        ) = self._parse_block(
            buffer,
            data,
            dataStart,
        )

        return self._finish_alignment(
            records,
            strands,
            annotations,
            score,
            parser,
        )