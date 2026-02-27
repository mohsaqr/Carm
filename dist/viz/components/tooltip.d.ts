/**
 * Hover tooltip component.
 * Shows a floating tooltip near the cursor with formatted content.
 */
import type { CarmTheme } from '../themes/default.js';
export declare function showTooltip(event: MouseEvent, content: string, theme?: CarmTheme): void;
export declare function hideTooltip(): void;
export declare function formatTooltipRow(label: string, value: string | number): string;
//# sourceMappingURL=tooltip.d.ts.map