import { slides } from "./deck-data.mjs";
import { addConfiguredSlide } from "./shared.mjs";

export default async function slide20(presentation, ctx) {
  return addConfiguredSlide(presentation, ctx, slides[19]);
}
