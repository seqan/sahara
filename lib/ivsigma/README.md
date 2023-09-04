# IVSigma
A c++20 library for bioinformatics. It provides functionality to work on nucleotide alphabets (e.g. dna4, dna5, rna4, ...), amino acid based alphabets and many more.

**IVSigma** stands for 4th version, and **Sigma** as the Greek letter Σ which is typically used to represent
alphabets as used in **bioinformatics**, such as Σ={A, C, G, T, N}

## Documentation/API Reference
Checkout the [documentation](https://iv-project.github.io/IVSigma/) for more information.

## Example

Here a quick example on how to convert dna5 data into rank space:

```c++
int main() {
    auto input  = std::string{"ACGnACGt"};
    auto output = ivs::convert_char_to_rank<ivs::dna5>(input);
    fmt::print("{} => {}\n", input, output);
    // prints: "ACGnACGt => [0, 1, 2, 4, 0, 1, 2, 3]"
}
```
