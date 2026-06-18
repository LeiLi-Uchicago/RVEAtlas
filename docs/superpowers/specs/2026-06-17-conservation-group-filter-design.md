# Conservation Group Filter Design

## Goal

Reorganize the Conservation tab controls so the gene selector is visually primary and optional group filtering is clearly separated, disabled by default, and accompanied by a performance warning when enabled.

## User Interface

- Keep the existing descriptive help text in the control card.
- Move the `Gene` selector to the left side of the main control row.
- Add a bordered `Filter by group` box on the right side.
- Place an enable/disable switch in the filter box header.
- Place `Group by` and `Group` side by side beneath the switch.
- Use the approved compact side-by-side layout (visual option A).
- On initial app load, the switch is off and both selectors are disabled and visually muted.
- Keep the existing minimum-sequence, plot-font-size, download-format, and download controls below the main row.

## Interaction

- With filtering disabled, entropy calculations use all groups.
- Turning the switch on:
  - enables `Group by` and `Group`;
  - displays a modal dialog every time it is enabled;
  - warns: “Real-time calculations using a group filter may take longer. Please wait while the results are recalculated.”
- The modal has an `OK` button and does not block calculation after dismissal.
- Turning the switch off:
  - disables both filter selectors;
  - recalculates using all groups;
  - preserves the current selector values for reuse if filtering is enabled again.
- Existing entropy loading and recalculation feedback remains in place.

## Server Behavior

- Add a Boolean Conservation filter input with a default value of `FALSE`.
- Derive an effective group selection:
  - filter enabled: use the selected `ent_group`;
  - filter disabled: use `"All"`.
- Use the effective group selection for entropy data retrieval, plot title text, variable-site summaries, and other Conservation outputs that currently read `ent_group`.
- Continue updating available `Group by` and `Group` choices while the controls are disabled so they are ready when enabled.
- Avoid resetting the preserved group selection solely because filtering was disabled.

## Styling and Accessibility

- The switch label clearly identifies its purpose as `Filter by group`.
- Disabled selectors retain their labels and use native disabled behavior plus subdued box styling.
- The filter box remains readable at narrower widths by allowing the two fields to stack responsively.
- The modal title and message explain the performance impact without implying an error.

## Testing

- Verify the filter switch defaults to off.
- Verify both group selectors are disabled when the switch is off.
- Verify calculations use all groups while disabled, regardless of preserved selector values.
- Verify enabling the switch enables both selectors.
- Verify the warning modal appears each time the switch changes from off to on.
- Verify disabling the switch preserves the selected values.
- Verify enabled filtering uses the selected group in entropy calculations and labels.
- Verify the layout stacks cleanly at narrow viewport widths.
- Run syntax checks and the project's available automated tests, then visually verify the Conservation tab in the browser.

## Scope

This change is limited to the Conservation tab control layout and group-filter activation behavior. It does not alter entropy formulas, database queries, or controls on other tabs.
