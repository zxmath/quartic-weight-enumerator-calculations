# Quartic weight-enumerator calculations

SageMath 10.9 code for Xu Zhuang, *The Projective Reed–Muller Code from Quartic Curves in the Projective Plane*.

| Part | Notebook | Contents |
| --- | --- | --- |
| Reducible quartics | [reducible.ipynb](notebooks/reducible.ipynb) | Incidence counts and R_z(q). |
| Genus 0 | [genus_0.ipynb](notebooks/genus_0.ipynb) | Profile counts and C_d(q). |
| Genus 1 | [genus_1.ipynb](notebooks/genus_1.ipynb) | B_d(E), weighted trace counts, and class numbers. |
| Genus 2 | [genus_2.ipynb](notebooks/genus_2.ipynb) | H_q, V_q, W_q and hyperelliptic-model sums. |
| Genus 3 | [genus_3.ipynb](notebooks/genus_3.ipynb) | S_q(n) for q=5,7 from rational-point normal forms. |
| Weight enumerator | [weight_enumerator.ipynb](notebooks/weight_enumerator.ipynb) | Assemble the contributions using Theorem 4.1. |

## Running

Run the notebooks in table order with a SageMath 10.9 kernel. Use `notebooks` as the working directory. Calculations export JSON files to `results/`, which the weight-enumerator notebook reads.

From the project root:

```sh
for notebook in reducible genus_0 genus_1 genus_2 genus_3 weight_enumerator
do
  sage -python -m jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=1800 --ExecutePreprocessor.record_timing=False "notebooks/$notebook.ipynb" || exit 1
done
```

The genus-one calculation supports prime powers of characteristic greater than three. The genus-two model sum supports primes greater than three and evaluates q=5,7. It includes all squarefree degree-five and degree-six polynomials, with every nonzero leading coefficient. The genus-three census uses SageMath finite fields, polynomial arithmetic, and linear algebra. No separate compiler is needed. Sources and normalization factors are stated beside the calculations.

## Examples

The saved q=5,7 examples are complete weight enumerators. Their coefficients sum to q^15. The smooth counts use rational-point normal forms and the classified pointless quartics. The q=5 distribution agrees term by term with Sutherland's database.
