# Problem

Consider an NxM maze with partitions placed in between some maze cells. It is
easy to see that there are 2^(N*(M-1)+M*(N-1)) unique mazes, ignoring
rotational and mirror symmetries. How many of these mazes can be navigated
from the top left corner to the bottom right corner?

Sample 4x3 maze that can be navigated from the top left corner to the bottom
right corner:

```
  +---+---+---+---+
  | S---+     |   |
  +---+ | +---+   +
  |     +---+ |   |
  +   +---+ | +---+
  |       | +---F |
  +---+---+---+---+
```

Sample 4x3 maze that cannot be navigated from the top left corner to the bottom
right corner:

```
  +---+---+---+---+
  | S             |
  +---+   +---+   +
  |       |   |   |
  +   +---+   +---+
  |   |         F |
  +---+---+---+---+
```

# Solution

The approach is to generate all possible paths from the top left corner to the
bottom right and count all partition placements that do not interfere with
those paths.  Each maze can be represented as a N*(M-1)+M*(N-1) bit mask with
1s in place of partitions. Each path through the maze forces some of the bits
in that maze bit mask to remain zeros while all the remaining bits can take any
value.

For example, assume that partitions are numbered as follows in the bit mask:

```
  +---+---+---+---+
  |   0   1   2   |
  +-9-+10-+11-+12-+
  |   3   4   5   |
  +13-+14-+15-+16-+
  |   6   7   8   |
  +---+---+---+---+
```

Consider this particular path:

```
  +---+---+---+---+
  | S---+         |
  +   + | +   +   +
  |     +---+     |
  +   +   + | +   +
  |         +---F |
  +---+---+---+---+
```

This path goes through partitions 0, 10, 4, 15 and 8. Therefore, the mask for
this path becomes 10111101011101110. There are 12 one bits in this mask, so
there are 2^12=4096 possible mazes that allow for this path through the maze.

Unfortunately we cannot simply generate all possible paths through the maze and
add up all numbers. That will count certain mazes several (potentially many)
times. As an example, an empty maze works for any path.

One approach is to simply generate all paths, convert them into masks, then run
through all possible mazes and test each maze to see if it satisfies at least
one path. This can be implemented in C using code similar to this:

```
  for (int i = 0; i < 1 << (n*(m-1)+m*(n-1)); i++) {
      if ((i & PATH_MASK_1) == 0 || (i & PATH_MASK_2) == 0 || ...) {
          counter++;
      }
  }
```

This becomes impractical when the size of the maze becomes larger than
approximately 4x4. For 5x5 maze you have to iterate over more than one trillion
(2^40) mazes and check each one against 8152 masks.

If all path masks were completely random and uniformly distributed, iterating
over all mazes and testing them would be the best way to go (see Boolean
Satisfiability problem). But they are not uniformly distributed, so there is a
better way.

For each path mask we can generate a stream of maze masks that satisfy this path
mask, then do a multiway merge between all those streams. We don't even have to
materialize all those streams. Since the stream can be produced already sorted,
we can use heap to perform the multiway merge as the streams are being
generated.
