import { GoogleGenAI } from "@google/genai";

const apiKey = process.env.GEMINI_API_KEY;
const model = process.env.GEMINI_MODEL || "gemini-3.5-flash";

if (!apiKey) {
  console.error(
    "GEMINI_API_KEY is not available. Set it as a Windows user environment variable, open a new terminal, and run this check again."
  );
  process.exit(2);
}

const client = new GoogleGenAI({
  apiKey,
  httpOptions: { apiVersion: "v1" }
});

try {
  const response = await client.interactions.create({
    model,
    input:
      "Reply with exactly: Gemini API connection successful. Do not add other text."
  });
  console.log(`Model: ${model}`);
  console.log(response.outputText || response.output_text || response.text);
} catch (error) {
  console.error("Gemini API connection failed.");
  console.error(error?.message || error);
  process.exit(1);
}
