# Pop-up dome generator

Inspired by Duncan Bimingham: <https://www.youtube.com/watch?v=2STS0POwB7g>

## Usage

The examples in the `examples` folder correspond to the following commands
(`'figsave', true` arguments omitted for brevity):

```matlab
% equivalent to dome(12, 14, 'earth', 'segmentdist', 'parallel')
dome()
dome(16, 18, {}, 'segmentdist', 'regular')
dome(8, 10, {}, 'segmentdist', 'maxvolume')
```

