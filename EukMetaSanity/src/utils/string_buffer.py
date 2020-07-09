from array import array


class StringBuffer:
    def __init__(self, output: str = None, buf_size: int = int(1e6), initial: str = ""):
        assert buf_size > 0
        assert isinstance(initial, str)
        assert output is None or isinstance(output, str)
        # Initialize buffer and position in buffer
        self._initialize_buffer(buf_size)
        self._buf_size = buf_size
        # Open output buffer
        self._output = None
        if output is not None:
            self._output = open(output, "w")
        # Add initial value, if present
        self._add_to_buffer(initial)

    def _initialize_buffer(self, buf_size: int):
        assert isinstance(buf_size, int)
        self._data = array("i", [0] * buf_size)
        self._pos = 0

    def _to_str(self) -> str:
        return "".join(chr(_char) for _char in self._data[:self._pos])

    @property
    def buffer(self) -> str:
        return self._to_str()

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
        self._initialize_buffer(self._buf_size)

    def __iadd__(self, other):
        self._add_to_buffer(other)
        return self

    def __repr__(self):
        return self._to_str()
