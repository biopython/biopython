typedef void (*Array_get_mapping_buffer_signature)(PyObject* self, Py_buffer* buffer);

typedef struct {
    PyObject* alphabet;
    Py_buffer mapping;
} Fields;
