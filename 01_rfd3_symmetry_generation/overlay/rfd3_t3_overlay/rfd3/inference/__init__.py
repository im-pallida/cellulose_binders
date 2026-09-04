# Sparse RFD3 overlay package.
# Files present here override original RFD3.
# Missing modules fall back to the installed RFD3 package.
try:
    from pkgutil import extend_path
    __path__ = extend_path(__path__, __name__)
except Exception:
    pass
