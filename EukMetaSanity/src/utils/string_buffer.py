from array import array


class StringBuffer:
    def __init__(self, output_path=None, buf_size=1e6, initial=""):
        assert buf_size > 0
        assert isinstance(initial, str)
        assert output_path is None or isinstance(output_path, str)
        # Initialize buffer and position in buffer
        self._buf_size = int(float(buf_size))
        self._initialize_buffer()
        # Open output buffer
        self._output = None
        if output_path is not None:
            self._output = open(output_path, "w")
        # Add initial value, if present
        self._add_to_buffer(initial)

    def _initialize_buffer(self):
        self._data = array("i", [0] * self._buf_size)
        self._pos = 0

    def _to_str(self) -> str:
        return "".join(chr(_char) for _char in self._data[:self._pos])

    @property
    def buffer(self) -> str:
        return self._to_str()

    @property
    def buf_size(self) -> int:
        return self._buf_size

    @buf_size.setter
    def buf_size(self, _size: int):
        self._flush_buffer()
        self.buf_size = _size
        self._initialize_buffer()

    def _add_to_buffer(self, _value: str):
        assert isinstance(_value, str)
        if self._pos + len(_value) >= self._buf_size:
            self._flush_buffer()

        for i, char in enumerate(_value):
            self._data[self._pos + i] = ord(char)
        self._pos += len(_value)

    def _flush_buffer(self):
        if self._output is not None:
            self._output.write(self._to_str())
        self._initialize_buffer()

    def __iadd__(self, other):
        self._add_to_buffer(other)
        return self

    def __repr__(self):
        return self._to_str()

    def add(self, value: str):
        self._add_to_buffer(value)
