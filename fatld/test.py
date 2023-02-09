from design import Design

if __name__ == "__main__":
    D = Design(32, 2, [5, 7, 11])
    print(D)
    print(D.array[0:5, :])
