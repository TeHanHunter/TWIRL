import { slides } from "./deck-data.mjs";
import { addConfiguredSlide } from "./shared.mjs";

export default async function slide23(presentation, ctx) {
  return addConfiguredSlide(presentation, ctx, slides[22]);
}
