# S56 Candidate Tensor Training Smoke

- label mode: `synthetic_label_smoke`
- synthetic label smoke: `True`
- training rows: `128`
- classes: `synthetic_bright, synthetic_faint`
- split counts: `{'test': 26, 'train': 76, 'validation': 26}`
- device: `cuda:0`
- torch: `2.11.0+cu128` / CUDA `12.8`
- final metrics: `{'epoch': 8.0, 'train_loss': 0.6935754292889645, 'train_accuracy': 0.5, 'train_n': 76.0, 'validation_accuracy': 0.5, 'validation_n': 26.0, 'test_accuracy': 0.5, 'test_n': 26.0}`

This is a bounded infrastructure smoke. Synthetic-label runs are not scientific models.
