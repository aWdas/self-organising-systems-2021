
def build_batch_list(prefix, start, end):
    print(" ".join([f"{prefix}{i}" for i in range(start, end)]))

build_batch_list("AF3920", 63, 83)