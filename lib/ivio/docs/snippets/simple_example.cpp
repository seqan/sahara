#include <ivio/ivio.h>

int main() {
    auto reader = ivio::fasta::reader{{"somedata.fasta"}};
    for (auto record_view : reader) {
        std::cout << "id: " << record_view.id << "\n";
        std::cout << "seq: " << record_view.seq << "\n";
    }
}
