# New Screenshots Directory

Place new screenshots and figures here for your updated book chapters.

### Recommended Usage in Quarto (`.qmd`) Chapters

To insert a new screenshot into any chapter:

```markdown
![Caption describing the figure or app interface](images/new_screenshots/your_screenshot.png){#fig-your-id width="85%" fig-align="center"}
```

You can then cross-reference it in text using `@fig-your-id`.
