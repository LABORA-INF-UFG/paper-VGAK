links = [(0, 1), (0, 2), (0, 5), (0, 10), (0, 13), (0, 14), (1, 2), (1, 4), (1, 5), (2, 10), (3, 4), (3, 5), (3, 7), (3, 8), (3, 13), (4, 5), (5, 13), (6, 7), (6, 8), (6, 16), (7, 8), (8, 12), (8, 13), (8, 16), (9, 10), (9, 11), (9, 14), (9, 19), (9, 20), (10, 11), (10, 14), (11, 20), (12, 13), (12, 14), (12, 16), (12, 17), (12, 19), (13, 14), (14, 19), (15, 16), (15, 17), (15, 18), (16, 17), (17, 19), (18, 17), (18, 19), (18, 20), (19, 20)]


new_links = []
print(links)

for l in links:
    if l[0] == 0:
        new_l0 = 1
    if l[0] == 1:
        new_l0 = 9
    if l[0] == 2:
        new_l0 = 10
    if l[0] == 3:
        new_l0 = 12
    if l[0] == 4:
        new_l0 = 13
    if l[0] == 5:
        new_l0 = 4
    if l[0] == 6:
        new_l0 = 18
    if l[0] == 7:
        new_l0 = 19
    if l[0] == 8:
        new_l0 = 3
    if l[0] == 9:
        new_l0 = 5
    if l[0] == 10:
        new_l0 = 11
    if l[0] == 11:
        new_l0 = 15
    if l[0] == 12:
        new_l0 = 7
    if l[0] == 13:
        new_l0 = 14
    if l[0] == 14:
        new_l0 = 2
    if l[0] == 15:
        new_l0 = 20
    if l[0] == 16:
        new_l0 = 8
    if l[0] == 17:
        new_l0 = 21
    if l[0] == 18:
        new_l0 = 6
    if l[0] == 19:
        new_l0 = 16
    if l[0] == 20:
        new_l0 = 17

    if l[1] == 0:
        new_l1 = 1
    if l[1] == 1:
        new_l1 = 9
    if l[1] == 2:
        new_l1 = 10
    if l[1] == 3:
        new_l1 = 12
    if l[1] == 4:
        new_l1 = 13
    if l[1] == 5:
        new_l1 = 4
    if l[1] == 6:
        new_l1 = 18
    if l[1] == 7:
        new_l1 = 19
    if l[1] == 8:
        new_l1 = 3
    if l[1] == 9:
        new_l1 = 5
    if l[1] == 10:
        new_l1 = 11
    if l[1] == 11:
        new_l1 = 15
    if l[1] == 12:
        new_l1 = 7
    if l[1] == 13:
        new_l1 = 14
    if l[1] == 14:
        new_l1 = 2
    if l[1] == 15:
        new_l1 = 20
    if l[1] == 16:
        new_l1 = 8
    if l[1] == 17:
        new_l1 = 21
    if l[1] == 18:
        new_l1 = 6
    if l[1] == 19:
        new_l1 = 16
    if l[1] == 20:
        new_l1 = 17

    new_links.append((new_l0, new_l1))

print(len(links))
print(len(new_links))

print(new_links)
print(links)

A = [[0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
     [1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
     [0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
     [0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1],
     [1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
     [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0],
     [1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
     [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0],
     [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0]]