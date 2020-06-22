#ifdef _TORCH_ENABLED
#include <torch/script.h> // One-stop header.
#endif //_TORCH_ENABLED

#include <iostream>
#include <memory>

int main(int argc, const char* argv[]) {
#ifdef _TORCH_ENABLED
	if (argc != 2) {
		std::cerr << "usage: example-app <path-to-exported-script-module>\n";
		return -1;
	}

	torch::jit::script::Module module;
	try {
		// Deserialize the ScriptModule from a file using torch::jit::load().
		module = torch::jit::load(argv[1]);
	}
	catch (const c10::Error& e) {
		std::cerr << "error loading the model\n";
		return -1;
	}

	// Create a vector of inputs.
	std::vector<torch::jit::IValue> inputs;
	inputs.push_back(torch::ones({1, 95}) * 0.1);

	// Execute the model and turn its output into a tensor.
	at::Tensor output = module.forward(inputs).toTensor();
	std::cout << output[0][0].item<float>() << '\n';

	std::cout << "ok\n";
#endif //_TORCH_ENABLED
	return 0;
}