from abc import ABC, abstractmethod
from typing import Optional
import os

# Try to load .env file with robust error handling
# Skip if file system operations fail (e.g., OneDrive sync issues)
try:
    from dotenv import load_dotenv

    # Use current working directory instead of __file__ to avoid path issues
    import os as os_module

    cwd = os_module.getcwd()
    env_file = os_module.path.join(cwd, ".env")
    if os_module.path.exists(env_file):
        load_dotenv(dotenv_path=env_file)
    else:
        # Fallback: try loading from current directory (no path specified)
        load_dotenv()
except (OSError, Exception):
    # Silently fail - environment variables might be set elsewhere
    # This allows the code to work even if .env file can't be read
    pass


class LLMProvider(ABC):
    """Abstract base class for LLM providers."""

    @abstractmethod
    def generate_response(self, prompt: str, system_prompt: Optional[str] = None) -> str:
        """
        Generate a response from the LLM.

        Args:
            prompt: User's question/prompt
            system_prompt: Optional system message to set context

        Returns:
            LLM response as string
        """
        pass


class OpenAIProvider(LLMProvider):
    """OpenAI API provider implementation."""

    def __init__(self, api_key: Optional[str] = None, model: str = "gpt-3.5-turbo"):
        try:
            from openai import OpenAI
        except ImportError:
            raise ImportError("OpenAI package not installed. Run: pip install openai")

        self.api_key = api_key or os.getenv("OPENAI_API_KEY")
        if not self.api_key:
            raise ValueError(
                "OPENAI_API_KEY not found in environment variables. Set it in .env file or environment."
            )

        try:
            self.client = OpenAI(api_key=self.api_key)
        except Exception as e:
            raise ValueError(f"Failed to initialize OpenAI client: {e}")

        self.model = model

    def generate_response(self, prompt: str, system_prompt: Optional[str] = None) -> str:
        messages = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": prompt})

        try:
            response = self.client.chat.completions.create(
                model=self.model, messages=messages, temperature=0.7, max_tokens=500
            )
            return response.choices[0].message.content
        except Exception as e:
            raise RuntimeError(f"OpenAI API error: {e}")


class AnthropicProvider(LLMProvider):
    """Anthropic (Claude) API provider implementation."""

    def __init__(self, api_key: Optional[str] = None, model: str = "claude-3-haiku-20240307"):
        try:
            from anthropic import Anthropic
        except ImportError:
            raise ImportError("Anthropic package not installed. Run: pip install anthropic")

        self.api_key = api_key or os.getenv("ANTHROPIC_API_KEY")
        if not self.api_key:
            raise ValueError(
                "ANTHROPIC_API_KEY not found in environment variables. Set it in .env file or environment."
            )

        try:
            self.client = Anthropic(api_key=self.api_key)
        except Exception as e:
            raise ValueError(f"Failed to initialize Anthropic client: {e}")

        self.model = model

    def generate_response(self, prompt: str, system_prompt: Optional[str] = None) -> str:
        messages = [{"role": "user", "content": prompt}]

        try:
            response = self.client.messages.create(
                model=self.model,
                max_tokens=500,
                system=system_prompt or "You are a helpful assistant.",
                messages=messages,
            )
            return response.content[0].text
        except Exception as e:
            raise RuntimeError(f"Anthropic API error: {e}")


class HuggingFaceProvider(LLMProvider):
    """HuggingFace provider implementation (for local models or API)."""

    def __init__(
        self, api_key: Optional[str] = None, model: str = "mistralai/Mistral-7B-Instruct-v0.2"
    ):
        try:
            from huggingface_hub import InferenceClient
        except ImportError:
            raise ImportError(
                "HuggingFace Hub package not installed. Run: pip install huggingface-hub"
            )

        self.api_key = api_key or os.getenv("HUGGINGFACE_API_KEY")
        self.model = model

        try:
            if self.api_key:
                self.client = InferenceClient(token=self.api_key)
            else:
                # Try without API key (for some public models)
                self.client = InferenceClient()
        except Exception as e:
            raise ValueError(f"Failed to initialize HuggingFace client: {e}")

    def generate_response(self, prompt: str, system_prompt: Optional[str] = None) -> str:
        full_prompt = f"{system_prompt}\n\n{prompt}" if system_prompt else prompt

        try:
            response = self.client.text_generation(
                prompt=full_prompt, max_new_tokens=500, temperature=0.7
            )
            return response
        except Exception as e:
            raise RuntimeError(f"HuggingFace API error: {e}")


def get_llm_provider(provider_name: Optional[str] = None) -> LLMProvider:
    """
    Factory function to get the appropriate LLM provider.

    Args:
        provider_name: Name of provider ("openai", "anthropic", "huggingface")
                      If None, uses DEFAULT_LLM_PROVIDER from env

    Returns:
        LLMProvider instance

    Raises:
        ValueError: If provider is unknown or API key is missing
        ImportError: If required package is not installed
    """
    try:
        provider_name = provider_name or os.getenv("DEFAULT_LLM_PROVIDER", "openai").lower()
    except Exception as e:
        raise ValueError(f"Failed to read DEFAULT_LLM_PROVIDER: {e}")

    if provider_name == "openai":
        return OpenAIProvider()
    elif provider_name == "anthropic":
        return AnthropicProvider()
    elif provider_name == "huggingface":
        return HuggingFaceProvider()
    else:
        raise ValueError(
            f"Unknown provider: {provider_name}. Choose: openai, anthropic, huggingface"
        )
