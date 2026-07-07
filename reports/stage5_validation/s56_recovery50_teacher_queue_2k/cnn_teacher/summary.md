# S56 Recovery50 CNN Teacher

- classes: `planet_like, instrumental_or_systematic, stellar_variability`
- class counts: `{'instrumental_or_systematic': 1511, 'planet_like': 347, 'stellar_variability': 122}`
- split counts: `{'test': 392, 'train': 1203, 'validation': 385}`
- group split: `tic`
- device: `cuda:0`
- best profile: `cnn_shape_plus_bls`

| Profile | Metadata | Best Val Balanced Accuracy | Test Balanced Accuracy | Real Balanced Accuracy | Injected Balanced Accuracy |
|---|---:|---:|---:|---:|---:|
| `cnn_shape_only` | `False` | 0.854 | 0.836 | 0.690 | 0.457 |
| `cnn_shape_plus_bls` | `True` | 0.916 | 0.885 | 0.746 | 0.787 |
